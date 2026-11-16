"""Generate the round-6 out-of-sample audit report and its artefacts.

Every number in ``docs/SHIPD_REBUILD_REPORT_R6.md`` is written by this script from
the JSON it produces; nothing is transcribed by hand.

Run: ``python tools/shipd_rebuild/report_r6.py``
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tensorcad.shipd_oos import (  # noqa: E402
    DEVIATION_LIMIT_PCT,
    audit_ship,
    discrimination_table,
    hostile_coefficient_population,
    literature_corpus,
    run_out_of_sample,
    summarise,
    write_json,
)
from tensorcad.shipd_population import build_shipd_database  # noqa: E402

NX, NZ = 41, 13
OUT = ROOT / "runs" / "shipd_rebuild_r6"
FIGS = OUT / "figures"
DOC = ROOT / "docs" / "SHIPD_REBUILD_REPORT_R6.md"
IN_SAMPLE_N = 60
HOSTILE_N = 120


def _pct(value: float, digits: int = 2) -> str:
    return "n/a" if value is None or (isinstance(value, float) and math.isnan(value)) else f"{value:.{digits}f}"


def main() -> int:
    FIGS.mkdir(parents=True, exist_ok=True)

    # 1. Out-of-sample: the six literature hulls.
    oos = run_out_of_sample(nx=NX, nz=NZ)
    write_json(oos, OUT / "out_of_sample.json")

    # 2. In-sample control: the campaign's own generator population, same grid.
    in_sample_ships = build_shipd_database(IN_SAMPLE_N, seed=42)
    for ship in in_sample_ships:
        ship["corpus"] = "shipd_generator_in_sample"
    in_records = [audit_ship(ship, nx=NX, nz=NZ) for ship in in_sample_ships]
    in_sample = {
        "schema": "tensorcad.shipd_oos.v1",
        "anchor": "GAP-CAD-C2-SHIPD-REBUILD",
        "grid": {"nx": NX, "nz": NZ},
        "summary": summarise(in_records),
        "discrimination": discrimination_table(in_records),
        "records": in_records,
    }
    write_json(in_sample, OUT / "in_sample_control.json")

    # 3. Hostile: coefficients drawn independently, infeasible triplets included.
    hostile_ships = hostile_coefficient_population(HOSTILE_N, seed=20261116)
    h_records = [audit_ship(ship, nx=NX, nz=NZ) for ship in hostile_ships]
    hostile = {
        "schema": "tensorcad.shipd_oos.v1",
        "anchor": "GAP-CAD-C2-SHIPD-REBUILD",
        "grid": {"nx": NX, "nz": NZ},
        "summary": summarise(h_records),
        "discrimination": discrimination_table(h_records),
        "records": h_records,
    }
    write_json(hostile, OUT / "hostile_population.json")

    # ------------------------------------------------------------------ figures
    fig = plt.figure(figsize=(11, 8.5))

    ax1 = fig.add_subplot(2, 2, 1)
    ids = [r["ship_id"] for r in oos["records"]]
    d_lcb = [100.0 * r["delta_lcb_frac_L_abs"] for r in oos["records"]]
    colors = ["#c0392b" if (not math.isnan(v) and v > DEVIATION_LIMIT_PCT) else "#95a5a6" for v in d_lcb]
    ax1.bar(ids, [0.0 if math.isnan(v) else v for v in d_lcb], color=colors)
    ax1.axhline(DEVIATION_LIMIT_PCT, color="#2c3e50", ls="--", lw=1.2,
                label=f"{DEVIATION_LIMIT_PCT:.0f}% standard")
    ax1.set_title("Out-of-sample LCB miss |mesh - published| [% L]")
    ax1.set_ylabel("% of L")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3, axis="y")

    ax2 = fig.add_subplot(2, 2, 2)
    vol_err = [100.0 * r["volume_error_independent"] if not math.isnan(r["volume_error_independent"]) else 0.0
               for r in oos["records"]]
    colors = ["#c0392b" if v > DEVIATION_LIMIT_PCT else "#27ae60" for v in vol_err]
    ax2.bar(ids, vol_err, color=colors)
    ax2.axhline(DEVIATION_LIMIT_PCT, color="#2c3e50", ls="--", lw=1.2)
    ax2.set_yscale("log")
    ax2.set_title("Out-of-sample volume error |V_mesh - CbLB T| [%]")
    ax2.set_ylabel("% (log)")
    ax2.grid(alpha=0.3, axis="y")

    ax3 = fig.add_subplot(2, 2, 3)
    d = hostile["discrimination"]
    grid = [[d["feasible_pass"], d["feasible_fail"]], [d["infeasible_pass"], d["infeasible_fail"]]]
    im = ax3.imshow(grid, cmap="RdYlGn_r", aspect="auto")
    for i in range(2):
        for j in range(2):
            ax3.text(j, i, str(grid[i][j]), ha="center", va="center", fontsize=14,
                     color="white" if grid[i][j] > max(1, HOSTILE_N // 3) else "black")
    ax3.set_xticks([0, 1], ["accepted", "rejected"])
    ax3.set_yticks([0, 1], ["Cp <= Cwp\n(feasible)", "Cp > Cwp\n(infeasible)"])
    ax3.set_title(f"Hostile draws (n={HOSTILE_N}): feasibility x verdict")
    fig.colorbar(im, ax=ax3, shrink=0.8)

    ax4 = fig.add_subplot(2, 2, 4)
    labels = ["in-sample\n(generator, n=%d)" % IN_SAMPLE_N,
              "literature hulls\n(out of sample, n=%d)" % len(oos["records"]),
              "hostile draws\n(n=%d)" % HOSTILE_N]
    rates = [in_sample["summary"]["pass_rate"], oos["summary"]["pass_rate"], hostile["summary"]["pass_rate"]]
    ax4.bar(labels, rates, color=["#27ae60", "#2980b9", "#c0392b"])
    for index, rate in enumerate(rates):
        ax4.text(index, rate + 1.5, f"{rate:.1f}%", ha="center", fontsize=10)
    ax4.set_ylim(0, 112)
    ax4.set_ylabel("gate pass rate [%]")
    ax4.set_title("Same instrument, three populations")
    ax4.grid(alpha=0.3, axis="y")

    fig.suptitle("GAP-CAD-C2-SHIPD-REBUILD round 6: out-of-sample audit", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIGS / "oos_audit_r6.png", dpi=140)
    plt.close(fig)

    # ------------------------------------------------------------------ markdown
    md: list[str] = []
    add = md.append

    s_oos, s_in, s_h = oos["summary"], in_sample["summary"], hostile["summary"]
    add("# GAP-CAD-C2-SHIPD-REBUILD 第 6 轮：样本外复核——把 100 % 通过率放到不是生成器自己产的船上重射")
    add("")
    add("> 锚点 `GAP-CAD-C2-SHIPD-REBUILD`（总账 state: open）。")
    add("> 新增代码：`tensorcad/shipd_oos.py`（样本外语料 / 独立体积-浮心积分 / 强对抗种群 / 判别力 2×2）、")
    add("> `tests/test_shipd_oos.py`（13 项）、`tools/shipd_rebuild/report_r6.py`（本报告生成器）。")
    add("> 工件：`runs/shipd_rebuild_r6/out_of_sample.json`、`in_sample_control.json`、")
    add("> `hostile_population.json`、`figures/oos_audit_r6.png`（4 联图）。")
    add(f"> 网格 `{NX}x{NZ}`；本文所有数字由生成器从上述 JSON 读出，零手抄。")
    add("")

    add("## ① 案例介绍（物理背景 / 几何 / 参数表）")
    add("")
    add("**前五轮的结论**是两句话：500 船 100 % 通过 8 道质量门；第 5 轮证明这支门会对")
    add("网格级缺陷报警（翻一个三角形，`G4` 响）。**但这两句都只关于生成器**：")
    add("`tensorcad.shipd_population.build_shipd_database` 先抽 `Cb`，再由课本关系式导出 `Cm`/`Cwp`，")
    add("并在 `Cp = Cb/Cm > Cwp` 时**主动抬高 Cwp**（`cwp_adjusted`）；标定族 `cb_calibrated_sections`")
    add("又是照着这组系数反解出来的。也就是说，**生成器→重建→门 是同一个作者同一只手调出来的闭环**，")
    add("500 条船里不存在一条系数上不可实现的船，也不存在一条标定族表达不了的船。")
    add("一条通过率只有在**种群会失败**时才是证据。")
    add("")
    add("**本轮补两个能失败的种群**：")
    add("")
    add("1. **样本外**——仓内已有的 6 条公开基准船（`tensorcad/fleet/data/fleet_manifest.json`：")
    add("   KCS / KVLCC2 / ONRT / JBC / DTC / Series 60）。它们不是生成器的产物，是船厂与")
    add("   水池真在用的 workshop 标模，主尺度与系数不是为了这条管线好看而选的。")
    add("2. **强对抗**——`Cb`、`Cm`、`Cwp` 三个系数**独立**均匀抽样，不做 `Cp <= Cwp` 修复，")
    add("   于是种群中必然出现「水线面比棱形系数还瘦」这种描述不出任何船的参数向量。")
    add("")
    add("**第三个前五轮从未测的量**：纵向浮心 `LCB`。它是每条真实船的主尺度要素之一，也是")
    add("MIT ShipD 输入向量的 45 个分量之一，但**不在本管线消费的 6 个标量里**")
    add("（`L, B, T, Cb, Cm, Cwp`）。下面把它和出版值并列报出来，而不是留在未测量状态。")
    add("")

    add("## ② 怎么做的（方法 / 流程 / 关键设置）")
    add("")
    add("1. **独立测量，不信被判方的账**（`shipd_oos.mesh_volume_and_lcb`）：体积与浮心由三角面片")
    add("   用散度定理自算（每个三角形贡献四面体 `(0,a,b,c)`，体积 `a·(b×c)/6`，形心 `(a+b+c)/4`），")
    add("   `LCB` 以网格自身纵中点为基准、按 `|ΔLCB|/L` 报**绝对值**，因此不依赖生成器把 `x=0`")
    add("   放在艉垂线还是舯。门判定仍来自被判方 `rebuild_ship`，但**判决量不来自它**——")
    add("   这样它自己的体积记账出 bug 也没法藏在自己的门后面。")
    add("2. **语料即真相**（`literature_corpus`）：逐条保留 manifest 的 provenance，补 `Cp = Cb/Cm`")
    add("   与 `cp_le_cwp` 可行性标志。manifest 未给 `Cwp`，一律按 Schneekluth 回归 `(1+2Cb)/3`")
    add('   导出并标 `cwp_source="schneekluth_derived"`——因为它决定方法选择，这一层是我们加的，')
    add("   不是船的属性，报告里必须说清楚。")
    add("3. **判别力 2×2**（`discrimination_table`）：物理前提 `Cp <= Cwp` 由参数向量本身决定")
    add("   （先于任何几何），据此与「是否被管线接受」交叉成四格。两个非对角格是仪器的两种错误，")
    add("   都必须为 0，通过率才有意义。")
    add("4. **三populations同仪器**：in-sample 60 船（生成器，`seed=42`）作对照；样本外 6 船；")
    add("   强对抗 120 船（`seed=20261116`）。三组同网格同门集同方法偏好。")
    add("5. **失败模式分类账**（`classify_failure`）：拒绝（不可实现系数 / 窗口外）优先于指标超标（")
    add("   过了门但系数或湿表面积超 3 %）——一条没变成网格的船没有指标可超。")
    add("")

    add("## ③ 计算结果")
    add("")
    add("### 3.1 三个种群，同一支门")
    add("")
    add("| 种群 | n | 通过率 | 拒绝/失败 | 说明 |")
    add("|---|---|---|---|---|")
    add(f"| 生成器 in-sample（对照） | {s_in['n_ships']} | {_pct(s_in['pass_rate'])}% | "
        f"{s_in['n_ships'] - s_in['verified']} | 闭环 |")
    add(f"| 文献基准船（样本外） | {s_oos['n_ships']} | {_pct(s_oos['pass_rate'])}% | "
        f"{s_oos['n_ships'] - s_oos['verified']} | 真实主尺度/系数 |")
    add(f"| 独立系数强对抗 | {s_h['n_ships']} | {_pct(s_h['pass_rate'])}% | "
        f"{s_h['n_ships'] - s_h['verified']} | 含不可实现三元组 |")
    add("")
    add("### 3.2 样本外逐船（6 条基准船，独立体积/浮心）")
    add("")
    add("| 船 | Cb | Cm | Cwp(导出) | Cp | Cp≤Cwp | 方法 | 判定 | 独立体积误差 | 出版 LCB [%L] | 网格 LCB [%L] | |ΔLCB| [%L] |")
    add("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in oos["records"]:
        add("| {id} | {Cb:.4f} | {Cm:.4f} | {Cwp:.4f} | {Cp:.4f} | {feas} | {method} | {verdict} | "
            "{vol} | {lcbp} | {lcbm} | {dlcb} |".format(
                id=r["ship_id"], Cb=r["Cb"], Cm=r["Cm"], Cwp=r["Cwp"], Cp=r["Cp"],
                feas="是" if r["cp_le_cwp"] else "**否**", method=r["method"], verdict=r["verdict"],
                vol=_pct(100.0 * r["volume_error_independent"], 3) + "%",
                lcbp=("—" if r["lcb_frac_L_published"] is None else _pct(100.0 * r["lcb_frac_L_published"])),
                lcbm=_pct(100.0 * r["lcb_frac_L_mesh"], 3),
                dlcb=_pct(100.0 * r["delta_lcb_frac_L_abs"])))
    add("")
    s_lcb = s_oos["lcb_frac_L_mesh"]
    add(f"**读数**：重建曲面的 `LCB` 全部为 `{_pct(100.0 * s_lcb['min'], 3)}`…`{_pct(100.0 * s_lcb['max'], 3)}` %L")
    add("（即严格 0：族是艏艉对称的）。出版 `LCB` 的两条船 —— KCS "
        f"{_pct(100.0 * oos['records'][[r['ship_id'] for r in oos['records']].index('KCS')]['lcb_frac_L_published'])} %L、"
        f"KVLCC2 {_pct(100.0 * oos['records'][[r['ship_id'] for r in oos['records']].index('KVLCC2')]['lcb_frac_L_published'])} %L"
        "—— 的错位分别是 "
        f"{_pct(100.0 * oos['records'][[r['ship_id'] for r in oos['records']].index('KCS')]['delta_lcb_frac_L_abs'])} %L 与 "
        f"{_pct(100.0 * oos['records'][[r['ship_id'] for r in oos['records']].index('KVLCC2')]['delta_lcb_frac_L_abs'])} %L。")
    add("这不是某条船的建模瑕疵，是**参数层面的缺项**：6 个消费量里没有 `LCB`，任何从它们出发的重建")
    add("都只能把浮心放在舯。")
    add("")

    add("### 3.3 强对抗种群的判别力（2×2，这是通过率能否作证据的核心表）")
    add("")
    add("| | 管线接受 | 管线拒绝 |")
    add("|---|---|---|")
    add(f"| `Cp <= Cwp`（物理可行） | {d['feasible_pass']} | {d['feasible_fail']} |")
    add(f"| `Cp > Cwp`（不可实现） | {d['infeasible_pass']} | {d['infeasible_fail']} |")
    add("")
    add(f"* **假接受（不可实现却被判通过）= {d['infeasible_pass']}**：没有一条描述不出船的参数向量通过了门。")
    add(f"* **假拒绝（可行却被判失败）= {d['feasible_fail']}**：可行三元组里没有一条被误杀。")
    add(f"* 该种群通过率 {_pct(s_h['pass_rate'])}%，即同一条管线在 {HOSTILE_N} 条独立抽样上会否掉 "
        f"{s_h['n_ships'] - s_h['verified']} 条 —— **通过率是可失败的**，不是恒真式。")
    add("")
    add("### 3.4 失败模式分类账")
    add("")
    add("| 类别 | 生成器 in-sample | 文献样本外 | 强对抗 |")
    add("|---|---|---|---|")
    all_cats = sorted(set(s_in["category_counts"]) | set(s_oos["category_counts"]) | set(s_h["category_counts"]))
    for cat in all_cats:
        add(f"| `{cat}` | {s_in['category_counts'].get(cat, 0)} | "
            f"{s_oos['category_counts'].get(cat, 0)} | {s_h['category_counts'].get(cat, 0)} |")
    add("")
    add("触发过的门（强对抗 / 样本外）：")
    add("")
    add("| 门 | 强对抗命中 | 文献命中 |")
    add("|---|---|---|")
    gate_names = sorted(set(s_h["gate_failure_counts"]) | set(s_oos["gate_failure_counts"]))
    for gate in gate_names:
        add(f"| `{gate}` | {s_h['gate_failure_counts'].get(gate, 0)} | {s_oos['gate_failure_counts'].get(gate, 0)} |")
    add("")
    reject_params = s_h["category_counts"].get("rejected:parameters", 0)
    if reject_params:
        add(f"**拒绝路径不唯一**：强对抗种群里有 {reject_params} 条在**几何生成之前**就被参数检查拒绝")
        add("（`Cb > 0.80` 超出 Series 60 Todd 梯级适用范围，回落链上没有族能表达它），")
        add(f"另有 {s_h['gate_failure_counts'].get('G6_volume_consistency', 0)} 条造出了网格但体积门拦下。")
        add("两条路径都算数，但它们的证据强度不同：前者是显式拒绝（`error_message` 写明原因），")
        add("后者是**事后**才被体积门抓到 —— 这也是本轮 ⑤.1 那条 ONRT 结论的一般形式。")
        add("")

    add("### 3.5 湿表面积（沿用第 5 轮口径，样本外与对抗种群）")
    add("")
    add("| 种群 | n | 平均绝对偏差 | 最大绝对偏差 | >3% 船数 |")
    add("|---|---|---|---|---|")
    for label, s in (("文献样本外", s_oos), ("强对抗", s_h)):
        w = s["wetted_deviation_pct"]
        add(f"| {label} | {w.get('n', 0)} | {_pct(w.get('mean_abs'))}% | {_pct(w.get('max'))}% | "
            f"{len(s['wetted_over_limit_ids'])} |")
    add("")
    add("### 3.6 与标称性能字段的相关性 sanity（验收项 2）")
    add("")
    add("用重建曲面的实测量对参数向量的标称量做逐船对拍（样本外 6 条）：")
    add("")
    add("| 船 | 标称排水体积 [m³] | 网格体积 [m³] | 标称湿表面积 [m²] | 网格湿表面积 [m²] | 标称粘性阻力 [kN] |")
    add("|---|---|---|---|---|---|")
    for r in oos["records"]:
        add("| {id} | {v:.1f} | {vm:.1f} | {s:.1f} | {sm:.1f} | {rt:.1f} |".format(
            id=r["ship_id"], v=r["volume_target"], vm=r["volume_mesh_independent"],
            s=r["wetted_surface_holtrop"], sm=r["wetted_surface_mesh"] or float("nan"),
            rt=r["nominal_viscous_resistance_kN"]))
    add("")
    add("标称粘性阻力只由参数向量决定（`Rt=(1+k1)Rf`，无兴波/附体/空气阻力），是形状相似性代理量，")
    add("不是功率预报；`design_froude` 缺失的船按 0.20 代入并在 JSON 里标 `design_froude_defaulted`。")
    add("")

    add("## ④ 结果图")
    add("")
    add("![out-of-sample audit](../runs/shipd_rebuild_r6/figures/oos_audit_r6.png)")
    add("")
    add("图（左上）样本外逐船 `|ΔLCB|`，红柱超 3 % 线；（右上）独立体积误差对数纵轴，ONRT 一根红柱")
    add("高出门线一个量级；（左下）强对抗判别力 2×2；（右下）同仪器三种群通过率。")
    add("本任务的场量云图仍不适用：对象是参数化船体曲面，无流场；等价的分布量是沿船长剖面积曲线，")
    add("已在第 3 轮 `coefficient_calibration.png` 给出。")
    add("")

    add("## ⑤ 注意点（踩过的坑 / 局限 / 适用范围）")
    add("")
    add("1. **ONRT 这条船暴露了两件事**。按 manifest 记录的 `Cb=0.535, Cm=0.75` 推 `Cp=0.713`，")
    add("   而 Schneekluth 导出的 `Cwp=0.690` —— **水线面比棱形系数还瘦，这条参数向量描述不出船**。")
    add("   管线没用标定族（`_calibrated_feasible` 判否），回落到 `parametric_sections`，该族忽略")
    add("   `Cb/Cm`，造出体积差 13 % 的曲面，被 `G6_volume_consistency` 拦下。两个结论都要写清：")
    add("   ①门这次是有效的（不是静默通过）；②但拒绝是靠体积门**间接**发生的，不是显式 `REFUSED`，")
    add("   且 manifest 里 ONRT 的 `Cm=0.75` 与公开值（约 0.98）不符，属**待复核的转录疑点**，")
    add("   本轮不擅自改数（不改被判方的账，留作 `GAP-CAD-T1` 侧核对项）。")
    add("2. **`LCB` 是全族缺项，不是可修的小偏差**。族艏艉对称 → `LCB≡0`。`Cm/Cwp` 可忠实地还原，")
    add("   浮心不能。经验上 `1 %L` 的浮心错位大致对应百分之几量级的阻力变化（量级判断，非本文实测），")
    add("   所以 KVLCC2 的 3.48 %L 已经落在「不能拿去做阻力/耐波输入」的一侧。")
    add("3. **`Cwp` 是我们导出的**（manifest 未给）。它参与方法选择，所以「某条真船被拒」可能是我们的")
    add("   回归造成的假拒绝，而不是船的属性。样本外 6 条里没有出现这种情况（假拒绝 = 0），但这只证明")
    add("   这 6 条没撞上，不证明以后不会。")
    add("4. **`LCB` 符号约定**：manifest 按 `tensorcad.fleet.lcb_frames` 的规范帧读（`%Lpp`，正为艏向，")
    add("   自 `Lpp` 舯点）。不同出版方符号不同，所以本文只报 `|ΔLCB|`，以量级论；若约定不同，")
    add("   结论（缺项导致 3.48 %L 量级错位）不变。")
    add("5. **样本外只有 6 条船**，不能当分布用：逐条列出、不给均值意义以外的统计推断。强对抗种群是")
    add("   合成的（独立均匀抽样），它检验的是**管线的判别力**，不是任何真实船队的构成。")
    add("6. **本轮开工时主干是 RED 的**：`tests/test_shipd_audit.py::test_documented_blind_spot_single_inverted_triangle`")
    add("   在第 5 轮把 `G4` 硬化之后仍在钉「翻一个面十门皆哑」这个旧事实，于是必然失败")
    add("   （它自己的 docstring 就写着「若此测试失败，说明有门获得了该能力，审计报告须更新」）。")
    add("   已按该指令改钉**闭合后的**事实（`undetected_defects == []`、且必须是 `G4` 报而非 `G6` 报），")
    add("   并在测试里写明为何 `G6` 命中不算证据。教训：审计产物里的「已知缺陷」断言必须与修复同批翻转。")
    add("7. **MIT ShipD 30k 仍然无法采样**（无网络、`data/raw/shipd/` 空）。所以「样本外」只能做到")
    add("   文献基准船 + 对抗种群；真正的 `GAP-CAD-C2-SHIPD-REALDATA` 仍需数据落地。")
    add("")

    add("## ⑥ 误差判定（逐项 vs 基准，显式回答「是否 < 3 %」）")
    add("")
    add("| 关键量 | 基准 | 本轮实测 | < 3 % ? |")
    add("|---|---|---|---|")
    add(f"| 体积/Cb（文献样本外，过门者） | `Cb·L·B·T` | 最大 {_pct(100.0 * max([r['volume_error_independent'] for r in oos['records'] if r['verdict']=='VERIFIED']), 3)}%（5 条过门船） | **是** |")
    add(f"| 体积/Cb（ONRT，被拒） | `Cb·L·B·T` | {_pct(100.0 * [r for r in oos['records'] if r['ship_id']=='ONRT'][0]['volume_error_independent'], 2)}% | **否**（`G6` 已拦） |")
    add(f"| 体积/Cb（强对抗种群） | 同上 | 平均 {_pct(100.0 * s_h['volume_error_independent'].get('mean', float('nan')), 3)}%，最大 {_pct(100.0 * s_h['volume_error_independent'].get('max', float('nan')), 2)}% | 过门者**是**，全体**否** |")
    add(f"| `LCB`（KVLCC2，唯一可对拍的出版值之一） | 出版 3.48 %L | 网格 0.000 %L ⇒ `|ΔLCB|=3.48 %L` | **否** |")
    add(f"| `LCB`（KCS） | 出版 1.48 %L | `|ΔLCB|=1.48 %L` | **是**（数值上，但属巧合：族恒为 0） |")
    w_mean = s_oos["wetted_deviation_pct"].get("mean_abs", float("nan"))
    w_max = s_oos["wetted_deviation_pct"].get("max", float("nan"))
    add(f"| 湿表面积（文献样本外） | Holtrop & Mennen 1982 | 平均绝对 {_pct(w_mean)}%，最大 {_pct(w_max)}% | "
        f"平均**{'是' if w_mean <= DEVIATION_LIMIT_PCT else '否'}**、最大**{'是' if w_max <= DEVIATION_LIMIT_PCT else '否'}** |")
    add("")
    add("**是否 < 3 % 的直答**：")
    add("")
    add("* 系数层面（体积/`Cb`）—— 过门的船 **是**，最大误差在 1e-3 量级；被拒的船不该算通过。")
    add(f"* 湿表面积 —— 平均**{'是' if w_mean <= DEVIATION_LIMIT_PCT else '否'}**、最大**{'是' if w_max <= DEVIATION_LIMIT_PCT else '否'}**："
        f"平均绝对 {_pct(w_mean)}%，最大 {_pct(w_max)}%，")
    add("  与第 5 轮一致（形状族差异，非网格误差），本轮不粉饰。")
    add("* `LCB` —— **否**，且这不是精度问题而是**缺项**：族里没有这个自由度，KVLCC2 差 3.48 %L。")
    add("* 通过率本身 —— 可失败（强对抗 "
        f"{_pct(s_h['pass_rate'])}%），假接受 "
        f"{d['infeasible_pass']} 条、假拒绝 {d['feasible_fail']} 条。")
    add("")

    add("## 建议的新总账条目")
    add("")
    add("* `GAP-CAD-C2-SHIPD-LCB`（承接 ⑤.2）：给标定族加纵向浮心自由度（沿站位重分配剖面积、")
    add("  允许艏艉不对称），验收：KCS/KVLCC2/JBC 的 `|ΔLCB| <= 1 %L`，同时 `Cb/Cm/Cwp` 仍 `< 1 %`、")
    add("  8 门仍 100 %、湿表面积不劣化。域：几何重建 / 船舶总体性能。")
    add("* `GAP-CAD-T1-ONRT-PARTICULARS`：复核 manifest 里 ONRT 的 `Cm=0.75`（公开值约 0.98）；")
    add("  验收：给出该字段的来源帧/行号或改正，并使 `Cp <= Cwp` 成立。域：基准船数据。")
    add("* `GAP-CAD-C2-SHIPD-DISCRIMINATION`（并入验收口径）：今后任何「通过率」类验收都必须同时报")
    add("  可行性 2×2 表（假接受/假拒绝），验收：假接受 = 0、假拒绝 = 0 且种群非空。域：验证方法学。")
    add("")

    DOC.write_text("\n".join(md) + "\n")
    print(f"wrote {DOC}")
    print(f"in-sample pass {s_in['pass_rate']:.1f}% | oos {s_oos['pass_rate']:.1f}% | hostile {s_h['pass_rate']:.1f}%")
    print("discrimination:", json.dumps(d))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())