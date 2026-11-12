import base64, json, pathlib

figs = []
caps = {
    "fig1_lift_models.png": "Left: module lift correlation (solid) vs classical Cv-dependent Savitsky statement (dashed), beta=20 deg. Right: their relative disagreement vs speed coefficient; no single Cv reconciles them.",
    "fig2_equilibrium_vs_assumed.png": "Left: wetted length ratio solved by the two-equation closure vs the constant lambda=3 the legacy solve echoed back. Right: running trim for the same displacement and beam.",
    "fig3_resistance_paths.png": "Total resistance from the two paths; the gap (+21% to +56% legacy) exceeds the +/-10% band the module claimed.",
    "fig4_moment_closure.png": "Centre of pressure vs wetted length ratio for several Cv; crossings with the LCG line are the moment-equilibrium solutions the legacy path could not produce.",
}
for name, cap in caps.items():
    p = pathlib.Path("artifacts/savitsky_audit") / name
    figs.append({"name": name, "caption": cap,
                 "png_base64": base64.b64encode(p.read_bytes()).decode()})

md = """## 1. 案例介绍

审计对象：`tensorcad.savitsky`（Savitsky 1964 滑行艇经验法，2026-11-12 由 m16869 合并，#359）。
物理背景：常斜升角棱柱形滑行艇在滑行状态下的升力、湿长比 λ、纵倾角 τ 与阻力预测。
参考艇（本审计全程使用）：Δ=15 kN，B=2.5 m，L=12 m，β=20°，LCG=3.0 m（艉垂线起沿龙骨）。
水：ρ=1025.9 kg/m³，ν=1.1892e-6 m²/s。有效性盒（模块自述）：Fr=1–4，β=10–30°，L/B=3–6。

问题：该模块合并时 28 个测试全绿、示例自证通过，但没有任何测试断言「报告出来的湿长比 λ 是解出来的」。

## 2. 怎么做的（方法/流程）

1. 独立重走执行路径（不复跑作者测试）：读 `running_trim` → `planing_prediction` → `resistance_components` 的耦合关系。
2. 探针矩阵（Δ=5–40 kN × V=6–15 m/s × B=2–3 m × β=10–30°，72 组）检查 λ 是否随输入变化。
3. 检查升力闭合（CL·½ρV²B² vs Δ）、封长（L_w > L）处理、有效性盒是否被求值、失败诊断是否可分辨。
4. 把两条升力相关式放在同一网格上对比（模块式 vs 经典含 Cv 式），并搜索是否存在单一 Cv 使二者一致。
5. 只对能证明的缺陷动手修复；依赖原始论文的部分标 WITHHELD，不臆造系数。
6. 全部结论写成可执行断言：`examples/savitsky_audit/run_example.py`（10 条 claim，任一失效即非零退出）+ 42 个测试。

## 3. 计算结果

### 3.1 阻塞缺陷 F1：所谓「均衡湿长比」从未被求解

| Δ [kN] | V [m/s] | B [m] | β [deg] | τ [deg] | 返回 λ |
|---:|---:|---:|---:|---:|---:|
| 8 | 8 | 2.5 | 20 | 3.20 | 3.000 |
| 12 | 10 | 2.0 | 10 | 2.52 | 3.000 |
| 15 | 12 | 2.5 | 20 | 2.07 | 3.000 |
| 20 | 15 | 3.0 | 25 | 0.90 | 3.000 |
| 26 | 9 | 2.2 | 30 | 6.62 | 3.000 |
| 40 | 18 | 3.0 | 15 | 2.75 | 3.000 |

六个物理上互不相关的工况返回同一个 λ=3.000（τ 从 0.90° 到 6.62°）。原因：一个方程两个未知量，λ 停留在初值。

### 3.2 新能力：两方程闭合后的解（升力=重量，压力中心=LCG）

| V [m/s] | Fr | Cv | τ 解 [deg] | λ 解 | R_T 解 [N] | τ 旧路径 [deg] | R_T 旧路径 [N] |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 0.74 | 1.62 | 3.057 | 2.124 | 1998 | 5.106 | 2904 |
| 12 | 1.11 | 2.42 | 2.501 | 1.804 | 2807 | 2.074 | 3825 |
| 16 | 1.47 | 3.23 | 1.905 | 1.708 | 3959 | 1.094 | 5853 |
| 20 | 1.84 | 4.04 | 1.515 | 1.667 | 5488 | 0.666 | 8571 |

升力残差 4.7e-15（相对），压力中心残差 1.2e-10 m；λ 随速度单调下降（2.124→1.667）；LCG 后移使纵倾上升（LCG=4.5 m 时 τ=1.40°，LCG=2.5 m 时 τ=2.96°）。

### 3.3 其他已修缺陷

- 封长钳制：L/B=2 工况原先把 λ 从 3.0 静默钳到 2.0 并破坏升力平衡；现在钳制写入 warnings 且在该 λ 下重解纵倾，升力平衡残差 2.4e-16。
- 有效性盒：Fr/β/trim/L-B 逐项求值并写入每个结果的 `warnings`；盒内工况返回空 warnings（测试双向钉住）。
- 诊断：原来的「Newton iteration failed to converge」实为窗口越界，现在直接报「required trim 47.25 deg is outside the implemented window」。
- 速度扫描：`speed_sweep(failures=[...])` 记录被跳过工况，`len(results)+len(failures)==len(speeds)` 有断言。

## 4. 结果图

图 1：两条升力相关式对比与分歧随 Cv 的变化。图 2：λ、τ 的「解出 vs 假定」对比。图 3：阻力后果。图 4：力矩闭合 l_cp(λ,Cv)=LCG 图解。均为程序化直出（无手抄）。

## 5. 注意点（坑与局限）

1. **升力相关式互不相容（F6，WITHHELD）**：模块式 `CL=0.00222 τ^0.9 λ^0.67 (β/1000)^-0.33` 不含速度系数 Cv，却被标注为「1964 论文式(1)」；经典式为 `CL_β=τ^1.1(0.0120λ^0.5+0.0055λ^2.5/Cv²)`、`CL_0=CL_β−0.0065βCL_β^0.6`。在 β=20°、τ=2–10°、λ=2–4 网格上模块式低 22%–78%，且最优单一 Cv 仍有 88% 最坏分歧——两者不可能同为同一篇论文的式(1)。本环境无原文，故未改模块公式，只把矛盾钉成测试（改动会让未验证系数冒充文献值）。
2. 阻力项 `R_P/Δ = tanτ + 0.008/√τ` 中的 `0.008/√τ` 无法溯源；它在 τ=1° 时占 tanτ 项的 46%，恰是滑行艇工作区，需原文核对。
3. 未做的检查：与 Savitsky 原论文、与模型试验数据（Sottorf / Kapryan & Boyd / Fridsma）、与 ship-performance-platform 独立实现的交叉比对——均需外部资料或跨仓访问。
4. 默认路径仍是旧的单方程路径（λ 假定），因为把默认切到 WITHHELD 相关式超出证据支持；新路径经 `lcg_m=` 显式开启并在结果里标注 provenance。

## 6. 误差判定（逐项 vs 基准）

| 关键量 | 基准 | 实测 | 是否 <3% |
|---|---|---|---|
| 升力闭合 CL·½ρV²B² vs Δ（闭合路径） | 重量 | 相对残差 4.7e-15 | 是（≪3%） |
| 压力中心 vs LCG（闭合路径） | LCG | 1.2e-10 m | 是（≪3%） |
| 封长钳制后升力闭合 | 重量 | 2.4e-16 | 是 |
| 闭式纵倾解 vs 升力方程 | 解析 | 0.0e+00 | 是 |
| 模块式 vs 经典式升力系数 | 经典式 | −22% … −78% | 否（远超 15% 不可接受线） |
| 阻力：闭合路径 vs 旧路径 | 旧路径 | 旧路径高 +21% … +56% | 否 |

诚实结论：本审计中与外部基准（实验/文献）可比对的量一个也没有（原文与试验数据不可得）。所有「通过」都是内部一致性与守恒类断言；两条路径之间 20%–80% 的差异说明该模块当前不能被引用为 ±10% 阻力预报器。F6 显式 WITHHELD，不粉饰。

工件：commit 39301d2；`docs/savitsky_audit_m16954.md`；`artifacts/savitsky_audit/{results.json,report.md,fig*.png}`；`tests/test_savitsky_audit.py`（42 测试）。"""

payload = {
    "title": "tensorcad.savitsky 对抗性审计：均衡湿长比从未求解的缺陷与两方程闭合补齐 (m16954)",
    "summary": "审计发现 savitsky 模块的纵倾求解从未解湿长比 lambda（恒为初值 3.0），湿面积/摩擦阻力/总阻力/有效功率因此都是常数的函数；封长钳制静默破坏升力平衡、有效性盒从不求值、失败诊断误导、速度扫描静默丢工况。新增 equilibrium_trim 两方程闭合解出 (tau, lambda)，残差 4.7e-15 / 1.2e-10 m；默认路径保持兼容但如实标注 provenance 与 warnings；模块升力式与经典含 Cv 式分歧 22%-78% 判为 WITHHELD 并钉成测试。",
    "scenario_key": "tensorcad-savitsky-audit",
    "mission_id": 16954,
    "tags": ["tensorcad", "savitsky", "planing", "audit", "verification", "withheld"],
    "markdown": md,
    "figures": figs,
}
pathlib.Path("/tmp/report.json").write_text(json.dumps(payload))
print("report json bytes:", len(json.dumps(payload)))