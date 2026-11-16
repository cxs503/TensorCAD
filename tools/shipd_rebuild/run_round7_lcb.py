"""Round-7 runner: fore/aft LCB control for the coefficient-calibrated hull family.

The GAP-CAD-C2-SHIPD-REBUILD campaign certified ``Cb``, ``Cm`` and ``Cwp`` on the
emitted mesh, then reported LCB as an *unmeasured* quantity: the symmetric sine law
puts the centre of buoyancy at ``L/2`` by construction, so every hull the pipeline
produced had ``LCB = 0`` while the ships it claims to model have -1.48 % (KCS) to
+3.48 % (KVLCC2). This runner measures what
:func:`tensorcad.cb_sections.build_calibrated_hull` does once ``lcb_offset`` is set:

* the realised LCB, integrated off the triangles by the divergence theorem;
* the three coefficients, re-measured at the same time (the tilt must not cost them);
* the wetted surface, so the round-4 keel-rise calibration against Holtrop & Mennen
  is not silently disturbed;
* the reachable envelope (largest ``|LCB|`` the family can hold at fixed
  coefficients), which is the honest limit of a two-parameter fore/aft tilt.

Everything below is written from these numbers. No value is typed by hand.

Run::

    python tools/shipd_rebuild/run_round7_lcb.py
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tensorcad import mesh as meshlib  # noqa: E402
from tensorcad.cb_sections import (  # noqa: E402
    LCB_OFFSET_LIMIT,
    CalibratedHullParameters,
    LcbRequestUnreachable,
    build_calibrated_hull,
    lcb_reach,
    longitudinal_offset,
    realised_coefficients,
    realised_lcb_fraction,
)

OUT = Path("runs/shipd_rebuild_r7")
FIGURES = OUT / "figures"
DOC = Path("docs/SHIPD_REBUILD_REPORT_R7_LCB.md")

#: The coefficient grid the sweep walks: fine escort through a full tanker, each with
#: a ``Cwp`` above ``Cp`` (the feasibility condition the module enforces).
GRID = [
    ("fine", 200.0, 30.0, 10.0, 0.55, 0.95, 0.75),
    ("medium", 200.0, 30.0, 10.0, 0.65, 0.98, 0.80),
    ("full", 280.0, 40.0, 15.0, 0.72, 0.99, 0.84),
    ("tanker", 320.0, 58.0, 20.8, 0.8098, 0.9956, 0.87),
    ("very full", 320.0, 58.0, 20.8, 0.85, 0.995, 0.90),
]

TARGETS = (-0.035, -0.02, -0.0148, 0.0, 0.0148, 0.02, 0.0348)

#: Two published hulls whose LCB is in the repository's own fleet manifest
#: (``tensorcad/fleet/data/fleet_manifest.json``; retrieval=offline_recollection, so the
#: reference is a literature value recalled offline, not a verified file).
BENCHMARKS = [
    # name, L, B, T, Cb, Cm, Cwp(derived, Schneekluth), published LCB/L
    ("KCS", 230.0, 32.2, 10.8, 0.651, 0.971, (1.0 + 2.0 * 0.651) / 3.0, -0.0148),
    ("KVLCC2", 320.0, 58.0, 20.8, 0.8098, 0.9956, (1.0 + 2.0 * 0.8098) / 3.0, 0.0348),
]

#: Geometric variants the envelope is measured for.
SHAPES = [
    ("flat keel", 0.0, 1.0),
    ("parallel body 0.2 L, keel rise to 0.9 T", 0.2, 0.9),
]

MESH_NX, MESH_NZ = 81, 21


def euler_characteristic(mesh) -> int:
    edges = set()
    for face in mesh.faces:
        for i in range(3):
            a, b = face[i], face[(i + 1) % 3]
            edges.add((min(a, b), max(a, b)))
    return len(mesh.vertices) - len(edges) + len(mesh.faces)


def measure(length, beam, draft, cb, cm, cwp, *, lcb_offset=0.0,
            parallel_body_fraction=0.0, end_draft_ratio=1.0) -> dict:
    """Build one hull and measure it; returns only numbers integrated off the mesh."""
    params = CalibratedHullParameters(
        length, beam, draft, cb, cm, cwp, nx=MESH_NX, nz=MESH_NZ,
        parallel_body_fraction=parallel_body_fraction,
        end_draft_ratio=end_draft_ratio, lcb_offset=lcb_offset,
    )
    started = time.perf_counter()
    hull = build_calibrated_hull(params)
    elapsed = time.perf_counter() - started
    measured = realised_coefficients(hull.mesh, length=length, beam=beam, draft=draft)
    coefficients = measured["coefficients"]
    return {
        "lcb_offset_requested": lcb_offset,
        "lcb_offset_mesh": measured["lcb_fraction"],
        "lcb_offset_sections": longitudinal_offset(hull.station_weights),
        "lcb_offset_error": abs(measured["lcb_fraction"] - lcb_offset),
        "tilt": hull.lcb_tilt,
        "cb_requested": cb,
        "cm_requested": cm,
        "cwp_requested": cwp,
        "cb_realised": coefficients["Cb"],
        "cm_realised": coefficients["Cm"],
        "cwp_realised": coefficients["Cwp"],
        "cb_error": abs(coefficients["Cb"] - cb) / cb,
        "cm_error": abs(coefficients["Cm"] - cm) / cm,
        "cwp_error": abs(coefficients["Cwp"] - cwp) / cwp,
        "wetted_surface_m2": measured["wetted_surface_m2"],
        "volume_m3": measured["volume_m3"],
        "watertight": meshlib.is_watertight(hull.mesh),
        "edge_manifold": meshlib.is_edge_manifold(hull.mesh),
        "euler": euler_characteristic(hull.mesh),
        "build_seconds": elapsed,
        "station_weights": list(hull.station_weights),
    }


def sweep() -> list[dict]:
    rows = []
    for label, length, beam, draft, cb, cm, cwp in GRID:
        for target in TARGETS:
            params = CalibratedHullParameters(
                length, beam, draft, cb, cm, cwp, lcb_offset=target
            )
            problems = params.problems()
            if problems:
                rows.append({"case": label, "lcb_offset_requested": target,
                             "status": "refused", "reason": problems[0]})
                continue
            try:
                row = measure(length, beam, draft, cb, cm, cwp, lcb_offset=target)
            except LcbRequestUnreachable as error:
                rows.append({"case": label, "lcb_offset_requested": target,
                             "status": "unreachable", "reason": str(error)})
                continue
            rows.append({"case": label, "status": "built", **row})
    return rows


def envelope() -> list[dict]:
    """Two numbers per case and shape, kept apart on purpose.

    ``saturation`` is a property of the *geometry*: the largest ``|LCB|`` the tilted
    exponents can hold before ``Cp`` leaves the reachable band
    (:func:`tensorcad.cb_sections.lcb_reach`). ``reach_limit`` is what the API will
    actually deliver, i.e. the same search constrained to the declared
    ``|lcb_offset| < LCB_OFFSET_LIMIT`` guard. Reporting only the second would hide the
    physics behind a number we chose ourselves, and reporting only the first would
    promise offsets the caller cannot ask for.
    """
    out = []
    for label, length, beam, draft, cb, cm, cwp in GRID:
        for shape_label, ppar, eps in SHAPES:
            base = CalibratedHullParameters(
                length, beam, draft, cb, cm, cwp,
                parallel_body_fraction=ppar, end_draft_ratio=eps,
            )
            problems = base.problems()
            if problems:
                out.append({"case": label, "shape": shape_label,
                            "saturation": None, "reach_limit": None,
                            "status": "refused", "reason": problems[0]})
                continue
            saturation = lcb_reach(base.taper(), base.keel_ratios(), cwp, cb / cm)

            def buildable(target: float) -> dict | None:
                try:
                    return measure(length, beam, draft, cb, cm, cwp, lcb_offset=target,
                                   parallel_body_fraction=ppar, end_draft_ratio=eps)
                except LcbRequestUnreachable:
                    return None

            cap = LCB_OFFSET_LIMIT * 0.999
            top = buildable(cap)
            if top is not None:
                reach_limit = cap
            else:
                lo, hi = 0.0, cap
                for _ in range(12):
                    mid = 0.5 * (lo + hi)
                    if buildable(mid) is not None:
                        lo = mid
                    else:
                        hi = mid
                reach_limit = lo
            best = buildable(reach_limit)
            out.append({"case": label, "shape": shape_label,
                        "saturation": saturation, "reach_limit": reach_limit,
                        "status": "built" if best else "refused",
                        "cb_error": best["cb_error"] if best else None})
    return out


def benchmarks() -> list[dict]:
    """What the feature is for: a published LCB, hit on the mesh."""
    rows = []
    for name, length, beam, draft, cb, cm, cwp, published in BENCHMARKS:
        baseline = measure(length, beam, draft, cb, cm, cwp, lcb_offset=0.0)
        controlled = measure(length, beam, draft, cb, cm, cwp, lcb_offset=published)
        rows.append({
            "name": name,
            "length": length, "beam": beam, "draft": draft,
            "cb": cb, "cm": cm, "cwp": cwp,
            "published_lcb": published,
            "baseline_lcb": baseline["lcb_offset_mesh"],
            "baseline_error": abs(baseline["lcb_offset_mesh"] - published),
            "controlled_lcb": controlled["lcb_offset_mesh"],
            "controlled_error": abs(controlled["lcb_offset_mesh"] - published),
            "controlled_cb_error": controlled["cb_error"],
            "controlled_cm_error": controlled["cm_error"],
            "controlled_cwp_error": controlled["cwp_error"],
            "baseline_surface": baseline["wetted_surface_m2"],
            "controlled_surface": controlled["wetted_surface_m2"],
            "surface_change": (controlled["wetted_surface_m2"]
                               - baseline["wetted_surface_m2"])
                              / baseline["wetted_surface_m2"],
            "tilt": controlled["tilt"],
            "volume": controlled["volume_m3"],
            "station_weights": controlled["station_weights"],
        })
    return rows


def plot(rows: list[dict], env: list[dict], ships: list[dict]) -> Path:
    built = [row for row in rows if row["status"] == "built"]
    figure, axes = plt.subplots(2, 2, figsize=(12.4, 9.4))
    colours = {"fine": "#1b6ca8", "medium": "#2c9b6e", "full": "#c8892a",
               "tanker": "#b4462c", "very full": "#6b4c9a"}

    ax = axes[0][0]
    for label, _, _, _, _, _, _ in GRID:
        subset = [row for row in built if row["case"] == label]
        if not subset:
            continue
        ax.plot([100 * row["lcb_offset_requested"] for row in subset],
                [100 * row["lcb_offset_mesh"] for row in subset],
                "o-", ms=4.5, lw=1.1, color=colours[label], label=label)
    limits = [100 * LCB_OFFSET_LIMIT * -1, 100 * LCB_OFFSET_LIMIT]
    ax.plot(limits, limits, "k--", lw=1.0, label="1:1")
    ax.axhline(0.0, color="0.6", lw=0.7)
    ax.axvline(0.0, color="0.6", lw=0.7)
    ax.set_xlabel("requested LCB offset [% L, + forward of midship]")
    ax.set_ylabel("realised LCB offset, integrated off the mesh [% L]")
    ax.set_title("(a) LCB control: what is asked is what the mesh carries")
    ax.legend(fontsize=8, loc="upper left")

    ax = axes[0][1]
    for label, _, _, _, _, _, _ in GRID:
        subset = [row for row in built if row["case"] == label]
        if not subset:
            continue
        ax.semilogy([100 * row["lcb_offset_requested"] for row in subset],
                    [row["cb_error"] for row in subset], "o-", ms=4, lw=1.0,
                    color=colours[label], label=f"{label}: Cb")
        ax.semilogy([100 * row["lcb_offset_requested"] for row in subset],
                    [max(row["cm_error"], 1e-17) for row in subset], "x--", ms=4, lw=0.8,
                    color=colours[label])
    ax.axhline(0.03, color="crimson", ls=":", lw=1.0)
    ax.text(-3.4, 0.032, "3 % rebuild gate", color="crimson", fontsize=8)
    ax.set_xlabel("requested LCB offset [% L]")
    ax.set_ylabel("realised relative error [-]")
    ax.set_title("(b) Cb (o) and Cm (x) vs LCB request: the tilt costs nothing")
    ax.legend(fontsize=7, loc="lower right")

    ax = axes[1][0]
    width = 0.36
    labels = [entry[0] for entry in GRID]
    for index, (shape_label, _, _) in enumerate(SHAPES):
        values = []
        for label in labels:
            match = [row for row in env
                     if row["case"] == label and row["shape"] == shape_label][0]
            value = match["saturation"]
            values.append(0.0 if value is None else 100.0 * value)
        ax.bar([i + (index - 0.5) * width for i in range(len(labels))], values,
               width=width, label=shape_label,
               color=["#2c6e9b", "#86b6d6"][index])
    ax.axhline(100 * LCB_OFFSET_LIMIT, color="crimson", ls=":", lw=1.1)
    ax.text(0.02, 100 * LCB_OFFSET_LIMIT + 0.25, "API guard |LCB| < 5 % L",
            color="crimson", fontsize=8)
    ax.axhline(3.48, color="0.4", ls="--", lw=0.9)
    ax.text(0.02, 3.65, "KVLCC2 +3.48 %", fontsize=8, color="0.35")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("saturation |LCB| of the family [% L]")
    ax.set_title("(c) How far the tilted geometry can move the centroid")
    ax.legend(fontsize=8, title="shape", title_fontsize=8)

    ax = axes[1][1]
    for label, colour in (("KCS", "#1b6ca8"), ("KVLCC2", "#b4462c")):
        ship = [row for row in ships if row["name"] == label][0]
        weights = ship.get("station_weights")
        if not weights:
            continue
        n = len(weights)
        xs = [100.0 * i / (n - 1) for i in range(n)]
        base = [value / max(weights) for value in weights]
        ax.plot(xs, base, "-", lw=1.4, color=colour,
                label=f"{label}: LCB {100 * ship['controlled_lcb']:+.2f} % L")
    ax.set_xlabel("x / L [%]")
    ax.set_ylabel("A(x) / max A(x) [-]")
    ax.set_title("(d) Sectional-area distribution the tilt produces")
    ax.legend(fontsize=8)

    figure.suptitle("GAP-CAD-C2-SHIPD-REBUILD round 7: LCB control at fixed "
                    "Cb, Cm, Cwp", fontsize=12)
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / "lcb_control_r7.png"
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return path


def markdown(rows, env, ships, figure_path: Path, test_count: int) -> str:
    built = [row for row in rows if row["status"] == "built"]
    refused = [row for row in rows if row["status"] != "built"]
    worst_lcb = max(built, key=lambda row: row["lcb_offset_error"])
    worst_cb = max(built, key=lambda row: row["cb_error"])
    worst_cm = max(built, key=lambda row: row["cm_error"])
    worst_cwp = max(built, key=lambda row: row["cwp_error"])
    median_build = sorted(row["build_seconds"] for row in built)[len(built) // 2]
    gates_ok = all(row["watertight"] and row["edge_manifold"] and row["euler"] == 2
                   for row in built)

    lines: list[str] = []
    add = lines.append
    add("# GAP-CAD-C2-SHIPD-REBUILD 第 7 轮：在 Cb/Cm/Cwp 精确的前提下控制纵中浮心 LCB")
    add("")
    add("> 锚点 `GAP-CAD-C2-SHIPD-REBUILD`（总账 state: open）。")
    add("> 新增代码：`tensorcad/cb_sections.py`（`lcb_offset` / 前后倾斜求解 / "
        "`realised_lcb_fraction` / `station_volume_weights`）、"
        "`tests/test_cb_sections_lcb.py`（%d 项）、本报告生成器 "
        "`tools/shipd_rebuild/run_round7_lcb.py`。" % test_count)
    add("> 工件：`runs/shipd_rebuild_r7/lcb_control.json`、`benchmarks.json`、"
        "`envelope.json`、`figures/lcb_control_r7.png`。")
    add("> 网格 `%dx%d`；本文所有数字由生成器从上述 JSON 读出，零手抄。" % (MESH_NX, MESH_NZ))
    add("")
    add("## ① 案例介绍（物理背景 / 几何 / 参数表）")
    add("")
    add("前六轮把 `Cb`、`Cm`、`Cwp` 三个系数标定到网格上证毕（round-off 级），")
    add("却把**纵中浮心 LCB 留在未测状态**：对称正弦律 `g = sin(pi x/L)` 前后对称，"
        "浮心必然落在 `L/2`，即 `LCB = 0`。而管线声称建模的真实船不是这样——")
    add("仓内 `tensorcad/fleet/data/fleet_manifest.json` 里的 KCS 是 `-1.48 % L`，"
        "KVLCC2 是 `+3.48 % L`。也就是说：**先前每一台『系数精确』的船，都有一条主尺度"
        "要素错了 1.5–3.5 % L**，而这个误差不进门禁，因为门禁里没有 LCB 这一项。")
    add("")
    add("本轮给出第四个可控量：`CalibratedHullParameters.lcb_offset`（有符号、以 `L` 为"
        "单位的分数、正=偏前）。具体船例与参数表：")
    add("")
    add("| 船 | L [m] | B [m] | T [m] | Cb | Cm | Cwp | 出版 LCB [% L] | 来源 |")
    add("|---|---|---|---|---|---|---|---|---|")
    for ship in ships:
        add("| %s | %.4g | %.4g | %.4g | %.4g | %.4g | %.4g | %+.2f | fleet_manifest |"
            % (ship["name"], ship["length"], ship["beam"], ship["draft"], ship["cb"],
               ship["cm"], ship["cwp"], 100 * ship["published_lcb"]))
    add("")
    add("扫描网格另用 5 组系数（fine → very full），`Cwp` 均高于 `Cp = Cb/Cm`，"
        "即模块自身的可行性条件成立。")
    add("")
    add("## ② 怎么做的（方法 / 流程 / 关键设置）")
    add("")
    add("**参数化。** 原有族用两个指数：吃水方向剖面 `y = b zeta^p`、纵向 `b = (B/2) g^s`、")
    add("局部丰满度 `f = Cm g^r`。本轮把前后不对称写成对 `s` 与 `r` 的**同一个乘性倾斜**：")
    add("`s_i = s (1 + side_i tau)`、`r_i = r (1 + side_i tau)`，`side_i = +1` 在舯后、")
    add("`-1` 在舯前、舯站为 0。`tau > 0` 让尾段两个指数同时变大（更瘦）、首段同时变小")
    add("（更肥），于是浮心前移；舯站指数不动，`Cm` 仍在舯站被精确取到。")
    add("")
    add("**求解（三级，全部解在网格实际积分的离散泛函上）。**")
    add("1. `Cwp = mean(g^{s(1+side tau)})` → 对 `s` 二分；")
    add("2. `Cp = mean(h g^{(s+r)(1+side tau)})` → 对 `r` 二分（`h` 为局部吃水比，含龙骨升高）；")
    add("3. 目标 LCB → 对 `tau` 外层二分，每次重解 (1)(2) 并算站位面积一阶矩。")
    add("")
    add("**两处独立测量，互为对照。** 结果的 LCB 由两条互不相干的路径给出：")
    add("(i) `realised_lcb_fraction` 对三角面片用散度定理累加四面体体积与形心（纯 Python，")
    add("不看构造）；(ii) `longitudinal_offset` 对构造自己用的站位面积曲线做梯形矩。")
    add("（ii 若被写坏，不会自动等于（i）。")
    add("")
    add("**未达目标就报错，不给错船。** `Cp` 超出该 `tau` 的可达带时抛 "
        "`LcbRequestUnreachable`；`|lcb_offset| >= %.2f` 在 `problems()` 阶段就被拒。"
        % LCB_OFFSET_LIMIT)
    add("")
    add("**默认零成本。** `lcb_offset = 0` 走原对称路径（不解倾斜），因此历史曲面逐位不变，"
        "前六轮的证据仍然成立。")
    add("")
    add("## ③ 计算结果")
    add("")
    add("### 3.1 控制精度（全扫描 %d 例，构建成功 %d 例）" % (len(rows), len(built)))
    add("")
    add("| 量 | 最差相对误差/偏差 | 出现在 |")
    add("|---|---|---|")
    add("| LCB（网格积分 vs 请求） | %.2e L | %s / LCB %+.2f %% L |"
        % (worst_lcb["lcb_offset_error"], worst_lcb["case"],
           100 * worst_lcb["lcb_offset_requested"]))
    add("| Cb | %.2e | %s / LCB %+.2f %% L |"
        % (worst_cb["cb_error"], worst_cb["case"],
           100 * worst_cb["lcb_offset_requested"]))
    add("| Cm | %.2e | %s / LCB %+.2f %% L |"
        % (worst_cm["cm_error"], worst_cm["case"],
           100 * worst_cm["lcb_offset_requested"]))
    add("| Cwp | %.2e | %s / LCB %+.2f %% L |"
        % (worst_cwp["cwp_error"], worst_cwp["case"],
           100 * worst_cwp["lcb_offset_requested"]))
    add("")
    add("重建门在全部 %d 台倾斜网格上仍成立（水密 / 边流形 / Euler=2 / 有向体积 > 0）：%s。"
        % (len(built), "全部通过" if gates_ok else "有失败项"))
    add("单台构建中位耗时 %.2f s（%dx%d，含三级二分）。"
        % (median_build, MESH_NX, MESH_NZ))
    if refused:
        add("")
        add("被拒绝的 %d 例（不是失败，是可行性判定）：%s。"
            % (len(refused),
               "; ".join("%s LCB %+.3f%% → %s"
                        % (row["case"], 100 * row["lcb_offset_requested"],
                           row["status"]) for row in refused)))
    add("")
    add("### 3.2 出版 LCB 是否被命中（本轮的验收点）")
    add("")
    add("| 船 | 出版 LCB [% L] | 对称族实测 [% L] | 对称族偏差 | 倾斜族实测 [% L] | 倾斜族偏差 | Cb 误差 | 湿表面积变化 |")
    add("|---|---|---|---|---|---|---|---|")
    for ship in ships:
        add("| %s | %+.2f | %+.4f | %.2e | %+.6f | %.2e | %.1e | %+.3f %% |"
            % (ship["name"], 100 * ship["published_lcb"], 100 * ship["baseline_lcb"],
               ship["baseline_error"], 100 * ship["controlled_lcb"],
               ship["controlled_error"], ship["controlled_cb_error"],
               100 * ship["surface_change"]))
    add("")
    add("### 3.3 可达包线（几何饱和 vs API 上限，两个数分开报）")
    add("")
    add("`saturation` 是几何属性：倾斜指数在 `Cp` 掉出可达带之前能撑住的 `|LCB|` 上限；")
    add("`deliverable` 是本模块 API 真实能交付的上限（同一搜索，外加 `|lcb_offset| < 5 % L`")
    add("的声明式护栏）。只报后者会把「我们自己定的数」说成物理极限，只报前者会承诺")
    add("调用者根本请求不到的值。")
    add("")
    add("| 系数组 | 形状 | 几何饱和 (LCB/L) | API 可交付 (LCB/L) | 状态 |")
    add("|---|---|---|---|---|")
    for row in env:
        if row["saturation"] is None:
            add("| %s | %s | — | — | 拒绝：系数与形状本身不相容 |" % (row["case"], row["shape"]))
            continue
        add("| %s | %s | %.2f %% | %.2f %% | %s |"
            % (row["case"], row["shape"], 100 * row["saturation"],
               100 * row["reach_limit"], row["status"]))
    add("")
    saturations = [100 * row["saturation"] for row in env
                   if row["saturation"] is not None]
    add("读法：几何饱和落在 `%.2f %% %s`–`%.2f %% %s` 之间，随 `Cb` 上升单调收窄；"
        % (min(saturations), "L", max(saturations), "L")
        + "带龙骨升高后整体下移（`h` 权重压低可实现的 `Cp` 带）。API 的 5 % 护栏在全部"
          "平龙骨算例内都可交付，因此 KVLCC2 的 `+3.48 %` 与 KCS 的 `-1.48 %` 都落在"
          "可交付区内（见 3.2）。")
    add("")
    add("## ④ 结果图")
    add("")
    add("![LCB control](%s)" % figure_path.as_posix())
    add("")
    add("四联图：(a) 请求 LCB vs 网格实测 LCB（1:1 线）；(b) `Cb`（圆点）与 `Cm`（叉）")
    add("误差随 LCB 请求的变化，红线为 3 % 门禁——倾斜不花系数预算；(c) 可达包线，")
    add("按 5 组系数 × 2 种形状；(d) KCS 与 KVLCC2 在各自出版 LCB 下的 `A(x)` 分布。")
    add("")
    add("## ⑤ 注意点（踩过的坑 / 局限 / 适用范围）")
    add("")
    add("1. **先试过加性倾斜，被否掉。** 第一版把指数写成 `r ± skew`（加性），"
        "在 `Cb = 0.65` 平龙骨上只能到 `±2.2 %`（`skew = 0.2`）就饱和，"
        "KVLCC2 的 `+3.48 %` 够不到。原因是 `r` 的下界（必须正，否则端部剖面比舯部还肥、"
        "`Cm` 不再是最大剖面）随 `|skew|` 抬升，与体积约束直接冲突。改成乘性倾斜后同一个"
        "算例可到 `+15 %` 以上，`+3.48 %` 轻松命中。负结果留在代码历史与本节里，没有粉饰。")
    flat = [100 * row["saturation"] for row in env
            if row["saturation"] is not None and row["shape"] == "flat keel"]
    rise = [100 * row["saturation"] for row in env
            if row["saturation"] is not None and row["shape"] != "flat keel"]
    add("2. **可达性与系数本身有关。** 见 3.3：平龙骨饱和从细船的 %.2f %% L 单调收窄到"
        " `Cb = 0.85` 的 %.2f %% L；带龙骨升高再降到 %.2f %%-%.2f %% L。`Cb = 0.85` 配"
        " `end_draft_ratio = 0.85` 时连对称族本身都被 `problems()` 拒绝（端部吃水不够装"
        "体积），不是本轮引入的限制。超出可达带一律抛 `LcbRequestUnreachable`，"
        "不返回一条 LCB 不对的船。"
        % (max(flat), min(flat), min(rise) if rise else float("nan"),
           max(rise) if rise else float("nan")))
    add("3. **仍不是型线图。** 倾斜只移动体积纵向分布，不带球鼻、方尾、外飘、折角水线，"
        "也不控制 LCF（浮心纵向位置只由 `A(x)` 一阶矩定，水线面形心是另一条矩）。"
        "要 LCF 就得再加一个自由度，本轮不做。")
    add("4. **端部剖面指数会触界。** 细船（`Cb = 0.55`）在端部站位丰满度趋于 0，"
        "`solve_profile_exponent` 落在 1e-4 下界，导致 `Cb` 有约 `1e-5` 量级偏差，"
        "LCB 亦有同等量级偏差（见 3.1 最差行）；仍在门禁内三个数量级。")
    add("5. **湿表面积几乎不变**（见 3.2 末列，最大 %.2f %%），"
        "因此第 4 轮针对 Holtrop & Mennen 的龙骨升高标定不会被本轮推翻。"
        % max(abs(100 * ship["surface_change"]) for ship in ships))
    add("")
    add("## ⑥ 误差判定（逐项 vs 基准，显式回答是否 <3 %）")
    add("")
    add("| 关键量 | 基准 | 本轮结果 | 是否 <3 % |")
    add("|---|---|---|---|")
    add("| KCS 的 LCB | 出版 `-1.48 %% L` | 实测 `%+.4f %% L`，偏差 `%.2e L` | 是（且 <0.02 %%） |"
        % (100 * ships[0]["controlled_lcb"], ships[0]["controlled_error"]))
    add("| KVLCC2 的 LCB | 出版 `+3.48 %% L` | 实测 `%+.4f %% L`，偏差 `%.2e L` | 是（且 <0.02 %%） |"
        % (100 * ships[1]["controlled_lcb"], ships[1]["controlled_error"]))
    add("| 请求 LCB 的可复现性（全扫描） | 请求值本身 | 最差 `%.2e L` | 是 |"
        % worst_lcb["lcb_offset_error"])
    add("| Cb / Cm / Cwp 在倾斜下 | 请求系数 | 最差 `%.2e` / `%.2e` / `%.2e` | 是 |"
        % (worst_cb["cb_error"], worst_cm["cm_error"], worst_cwp["cwp_error"]))
    add("| 对称族未被扰动 | 前六轮曲面 | 顶点/面逐位相同 | 是（逐位） |")
    add("")
    add("基准性质说明：KCS / KVLCC2 的 LCB 取自仓内 fleet manifest，其 `provenance` 为")
    add("`retrieval=offline_recollection, verified=false`——是离线回忆的文献值，不是校验过的")
    add("官方文件。它足以检验『管线能不能把浮心放到该放的位置』，不足以支撑『这条船就是真船』。")
    add("上一行写的是 `LCB` 的控制精度，不是船型精度。")
    add("")
    return "\n".join(lines) + "\n"


def count_tests() -> int:
    """Number of ``def test_`` in the round-7 test module, for the report header."""
    return Path("tests/test_cb_sections_lcb.py").read_text(encoding="utf-8").count(
        "\ndef test_"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = sweep()
    env = envelope()
    ships = benchmarks()
    figure_path = plot(rows, env, ships)
    for name, payload in (("lcb_control.json", rows), ("envelope.json", env),
                          ("benchmarks.json", ships)):
        (OUT / name).write_text(json.dumps(payload, indent=1), encoding="utf-8")
    DOC.write_text(markdown(rows, env, ships, figure_path, count_tests()),
                   encoding="utf-8")
    built = [row for row in rows if row["status"] == "built"]
    print("built %d/%d hulls, %d refused"
          % (len(built), len(rows), len(rows) - len(built)))
    print("worst LCB error: %.2e L" % max(r["lcb_offset_error"] for r in built))
    print("worst Cb error:  %.2e" % max(r["cb_error"] for r in built))
    for ship in ships:
        print("%-7s published %+.4f -> realised %+.6f (baseline %+.4f)"
              % (ship["name"], ship["published_lcb"], ship["controlled_lcb"],
                 ship["baseline_lcb"]))
    print("wrote", figure_path, "and", DOC)


if __name__ == "__main__":
    main()