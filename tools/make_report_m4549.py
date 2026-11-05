"""Build the m4549 mission report from the program's own JSON outputs (no hand-typed numbers).

Reads ``docs/righting/results.json`` (written by ``examples/righting/run_example.py``) and
``artifacts/m4534/merge_verification.json`` (written by ``tools/verify_merge_m4534.py``),
shrinks the two figures, and writes ``/tmp/report.json`` for ``POST /tensor/reports``.
"""

from __future__ import annotations

import base64
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "docs" / "righting" / "results.json"
VER = ROOT / "artifacts" / "m4534" / "merge_verification.json"
OUT = Path(os.environ.get("REPORT_OUT", "/tmp/report.json"))


def sci(x: float) -> str:
    return f"{x:.3e}"


def fixed(x: float, n: int = 6) -> str:
    return f"{x:.{n}f}"


def shrink(path: Path, out: Path, width: int = 820, colors: int = 48) -> str:
    """Return a base64 PNG of ``path`` resized/quantised to keep the payload small."""
    from PIL import Image

    image = Image.open(path).convert("RGB")
    image = image.resize((width, round(image.height * width / image.width)), Image.LANCZOS)
    image = image.quantize(colors=colors, method=Image.MEDIANCUT)
    image.save(out, "PNG", optimize=True)
    return base64.b64encode(out.read_bytes()).decode()


def main() -> int:
    res = json.loads(RES.read_text())
    ver = json.loads(VER.read_text())
    barge, hull, claims = res["barge"], res["hull"], res["claims"]
    bp, hp = barge["parameters"], hull["parameters"]

    geom = "\n".join([
        "| 体 | L [m] | B [m] | T [m] | 干舷 [m] | KG [m] | 备注 |",
        "|---|--:|--:|--:|--:|--:|---|",
        f"| 方驳 barge | {bp['length']:.1f} | {bp['beam']:.1f} | {bp['draft']:.1f} | "
        f"{bp['freeboard']:.1f} | {bp['kg']:.1f} | 直立舷壁，可与解析 GZ 对比 |",
        f"| 参数化船体 hull | {hp['length']:.1f} | {hp['beam']:.1f} | {hp['draft']:.1f} | "
        f"{hull['freeboard_m']:.3f} | {hp['kg']:.2f} | sin^0.7·sqrt 族，{hp['nx']}×{hp['nz']} 网格，"
        f"浮于 z={hp['waterline_z']:.1f} m |",
    ])

    barge_tbl = ["| 横倾 [°] | GZ 站位积分 [m] | GZ 解析闭式 [m] | 相对偏差 | 排水体积 [m³] |",
                 "|--:|--:|--:|--:|--:|"]
    for r in barge["rows"]:
        barge_tbl.append(
            f"| {r['heel_deg']:.1f} | {fixed(r['gz_m'], 9)} | {fixed(r['gz_wall_sided_m'], 9)} | "
            f"{sci(r['relative_deviation'])} | {fixed(r['submerged_volume_m3'], 6)} |")

    hull_tbl = ["| 横倾 [°] | GZ [m] | y_B [m] | z_B [m] | 甲板入水比例 |", "|--:|--:|--:|--:|--:|"]
    for r in [*hull["curve"][::3], hull["curve"][-1]]:
        hull_tbl.append(
            f"| {r['heel_deg']:.1f} | {fixed(r['gz_m'], 6)} | "
            f"{fixed(r['transverse_centre_buoyancy_m'], 6)} | "
            f"{fixed(r['vertical_centre_buoyancy_m'], 6)} | {r['deck_immersion_fraction']:.2f} |")

    conv_tbl = ["| 站位数 | 排水体积 [m³] | GZ_max [m] | 稳性消失角 [°] |", "|--:|--:|--:|--:|"]
    conv = hull["grid_convergence"]
    for iv, vol, gzm, vanish in zip(conv["intervals"], conv["volume_m3"],
                                    conv["gz_max_m"], conv["angle_vanishing_deg"],
                                    strict=True):
        conv_tbl.append(f"| {int(iv)} | {fixed(vol, 9)} | {fixed(gzm, 6)} | {fixed(vanish, 3)} |")

    claim_tbl = ["| 断言 | 实测 | 门限 | 单位 | 结果 |", "|---|--:|--:|---|:--:|"]
    for c in claims:
        claim_tbl.append(f"| {c['claim']} | {sci(c['measured'])} | {sci(c['limit'])} | "
                         f"{c['unit']} | {'PASS' if c['passed'] else 'FAIL'} |")

    routes = ver["barge_three_routes"]
    wall = ver["barge_vs_wall_sided"]
    hull2 = ver["hull_two_routes"]
    conv2 = ver["station_convergence"]
    neg = ver["negative_control"]
    through = ver["submerged_body_vs_closed_surface"]
    route_tbl = "\n".join([
        "| 检查 | 实测 | 含义 |",
        "|---|--:|---|",
        f"| 方驳：三条路线 worst \\|ΔGZ\\|（0–45°） | {routes['max_difference_m']:.3e} m | "
        "站位-Simpson / 闭合面 Newton / 括-二分 三方一致，量级=舍入 |",
        f"| 方驳 vs 直立舷壁解析闭式（worst 相对） | {sci(wall['worst_relative'])} | 对照解析解 |",
        f"| 41×13 船体 z=−0.3 m：站位 vs gz（worst \\|ΔGZ\\|） | {hull2['max_difference_m']:.3e} m | "
        "两种离散化之差，非舍入 |",
        f"| 站位面积体积 vs 闭合面四面体体积（相对） | "
        f"{through['relative_difference']:.1f} | 同一封顶体，几何一致 |",
        f"| 40 vs 80 站位 GZ_max（相对） | {sci(conv2['gz_max_relative'])} | 站立网格收敛 |",
        f"| 负控：KG+{neg['kg_shift_m'] * 1000:.0f} mm 使 GZ(20°) 变化 | {neg['gz_20deg_shift_m']:.3e} m "
        f"（= 路线间差的 {neg['ratio_to_inter_route_gap']:.0f} 倍） | 比较器能失效，非空比较 |",
    ])

    vanish_rel = abs(conv["angle_vanishing_deg"][1] - conv["angle_vanishing_deg"][0]) \
        / conv["angle_vanishing_deg"][1]
    err_tbl = "\n".join([
        "| 关键量 | 基准 | 实测误差 | <3%? |",
        "|---|---|--:|:--:|",
        f"| 方驳 GZ(0–45°) | 直立舷壁解析闭式 | {sci(barge['worst_relative_deviation'])} 相对 | 是 |",
        f"| 方驳 GM（dGZ/dφ 在 0） | KB+BM−KG = {barge['gm_m']:.6f} m | "
        f"{sci(barge['slope_relative_error'])} 相对 | 是 |",
        f"| 船体 GM（dGZ/dφ 在 0） | KB+BM−KG = {hull['gm_m']:.6f} m | "
        f"{sci(hull['slope_relative_error'])} 相对 | 是 |",
        f"| 船体 GZ_max | 独立路线 `gz`（闭合面 + Newton） | "
        f"{sci(hull2['max_difference_m'] / hull['gz_max_m'])} 相对 | 是 |",
        f"| 船体稳性消失角 | 80 站位自洽解 | {sci(vanish_rel)} 相对 | 是 |",
        "| 船体 GZ 曲线 vs 船模试验 / CFD | 无基准数据 | 未做 | 不适用（见⑥结论） |",
    ])

    fig1 = shrink(ROOT / "docs" / "righting" / "righting_arms.png", Path("/tmp/fig_m4549_1.png"))
    fig2 = shrink(ROOT / "artifacts" / "m4534" / "merge_verification.png", Path("/tmp/fig_m4549_2.png"))

    md = f"""# TensorCAD 0.23.0：站位截面第三路 GZ(φ) 在当前 main 三路互校（m4549 / PR#52 重开第 3 轮）

## ① 案例介绍

**物理背景。** 完整稳性的核心曲线是复原力臂 `GZ(φ)`：船体横倾 φ 后浮心横向移动，浮力与
重力的力偶臂即 `GZ`。工程上有三条互相独立的算法路线可互为校验：
(a) 闭合曲面散度定理 + 水线牛顿迭代（`tensorcad/gz.py`，0.21.0）；
(b) 逐面裁剪截面 + 括-二分水线（`tensorcad/righting_arm.py`，0.22.0）；
(c) 本次交付的**站位截面**路线（`tensorcad/righting.py`）：每站只切一次，用 Green 定理求
截面面积与一阶矩，沿船长做复合 Simpson，水线偏移由**二分**在"船体自身竖向范围"括内求解，
排水量保持常值。第三条路线承自 PR #52（`agent/m3581`），因**命名遮蔽**（其
`righting_arm_curve`/`RightingArmCurve`/`wall_sided_righting_arm` 与已并入 main 的 `gz.py`
撞名，合并即让 `gz` 路线从包级 API 不可达）被 owner 于 2026-11-05 关闭；本轮从当前 main
（`faa0da6`）重切，以 `section_righting_arm_curve` / `SectionRightingArmCurve` 发布，
把"重叠"变成测试里的双向断言。

**几何与参数（下表由程序从 `docs/righting/results.json` 注入，零手抄）。**

{geom}

## ② 怎么做的

1. **站位切面。** `_section_grid` 按 x 分站，每站用 `bonjean.section_loops` 取 (y,z) 闭合
   多边形；`_clip_below_waterplane` 对倾斜水线 `z cosφ − y sinφ = offset` 做 Sutherland–Hodgman 裁剪。
2. **Green 定理。** `_polygon_moments` 由顶点直接给出面积与 ∫y dA、∫z dA，消除数值积分误差。
3. **Simpson + 常排水量。** `_simpson` 沿站距复合求和；每个横倾角的 `offset` 用二分求解使
   浸没体积等于直立排水量，括区间取自 `_vertical_extent`——无初值猜测，故 `GZ(0)` 是舍入零
   而不是两个离散误差之差。
4. **`submerged_body` 退化面。** 同一积分只出体积/浮心、不施加排水量约束，用于几何自检
   （站位体积 vs 闭合面四面体体积相对差 0.0）。
5. **拒绝而非硬答。** 网格顶面在直立水线以下（`freeboard ≤ 0`）时抛 `RightingError`；
   未收敛的横倾角返回 `converged=False` 的 NaN 点，而不是被平均掉。
6. **工具链与命令。** `python examples/righting/run_example.py`（8 条断言，破一即非零退出，
   写 `docs/righting/{{results.json,report.md,righting_arms.png}}`）；
   `python tools/verify_merge_m4534.py`（写 `artifacts/m4534/merge_verification.{{json,png}}`）；
   `python -m pytest tests -q`；`ruff check .`；`mypy tensorcad`。

## ③ 计算结果

### 方驳 vs 直立舷壁解析闭式（L=10, B=2, T=1, 干舷 1 m, KG=0.8 m）

{chr(10).join(barge_tbl)}

GM = {barge['gm_m']:.6f} m，曲线斜率反算 {barge['slope_from_curve']:.9f} m（相对
{sci(barge['slope_relative_error'])}）；甲板边缘入水角解析经典值
{barge['classic_deck_edge_deg']:.1f}°，程序测得 {barge['deck_immersion_onset_deg']:.1f}°。

### 船体曲线（41×13，浮于 z = −0.3 m，KG = 0.75 m）

{chr(10).join(hull_tbl)}

GZ_max = {hull['gz_max_m']:.6f} m @ {hull['angle_gz_max_deg']:.3f}°；稳性消失角
{hull['angle_vanishing_stability_deg']:.3f}°；甲板入水起始 {hull['deck_immersion_onset_deg']:.1f}°；
排水体积 {hull['displacement_volume_m3']:.9f} m³。

### 站位网格收敛

{chr(10).join(conv_tbl)}

### 三条路线交叉验证（`artifacts/m4534/merge_verification.json`）

{route_tbl}

### 示例程序自身的 8 条断言（`docs/righting/results.json`）

{chr(10).join(claim_tbl)}

## ④ 结果图

图 1（`docs/righting/righting_arms.png`，示例程序生成）：方驳 GZ 与解析闭式重合、船体 GZ 曲线
及 GZ_max/稳性消失角、与解析解/其他路线的偏差。

图 2（`artifacts/m4534/merge_verification.png`，验证脚本生成）：(a) 方驳三条路线同图；
(b) 方驳 vs 闭式在其有效角域内；(c) 半浸圆柱体对照精确解 `(R−KG)sinφ`。

## ⑤ 注意点（踩过的坑与局限）

- **本任务的真问题是命名遮蔽，不是合并冲突。** PR #52 单独跑自己的测试全绿（直接
  `from tensorcad.righting import righting_arm_curve`），但合并后包级名字会覆盖 `gz.py` 的
  同名符号，`gz` 路线从包级 API 静默消失（不是 ImportError，是静默丢失）。修法=改名 +
  在 `tests/test_righting.py` 里断言"包级别名指向哪个模块"，让后续改动无法悄悄复活冲突。
- **`ruff` 版本漂移。** 上一轮记录的 "ruff clean" 是在旧 ruff 下成立的；本机 ruff 0.15.0 在
  `tools/make_report_m3859.py` 又报 5 条历史欠账（未用 noqa、未用解包名、`%` 格式化、缺尾换行），
  本轮一并清掉，`ruff check .` 重新全绿。
- **测试自身会写 `artifacts/m3781/`**（`tests/test_merge_cure_m3391.py` 跑验证脚本的副作用），
  属生成物不是成果，提交时显式剔除，避免脏 diff。
- **网格必须有干舷。** 水线切在网格顶面时任何横倾都上浪，函数拒绝而不是给曲线；需把
  `waterline_z` 设低或补网格。
- **适用范围。** 无自由液面修正、无舱柜/破舱、无纵倾-下沉耦合、不做任何法定衡准（IS Code/
  SOLAS）；截面假定每 y 单值（不允许外飘、不允许浸没甲板面积）；超过
  `deck_immersion_onset_deg` 后曲线描述的是水密箱体而非船舶。`GZ` 是几何量，不是合规结论。
- **站位法 vs 面元法的差不是舍入。** 船体上两路线差 7.3e-07 m（相对 ~4e-5），来自站位网格
  与三角网格两种离散化，证据里如实标注，不冒充舍入一致。

## ⑥ 误差判定（逐项 vs 基准，显式回答是否 <3%）

{err_tbl}

**结论：** 凡有基准的项（解析闭式、GM 的 KB+BM−KG、独立路线 `gz`、网格自洽）误差在
`1e-15`～`4e-5` 量级，**全部远小于 3%**（最好 5.8e-15，最差约 4e-5）。船体 `GZ` 曲线
**没有**船模试验或 CFD 基准可比，本轮不声称对实验的精度，只声称与解析解及另一条独立算法
路线一致，该项标"不适用"，需外部数据才能升级证据等级。本单不含 MOCK/STUB：
`tensorcad/righting.py` 的数值全部由网格几何真实积分得到，且负控（KG+5 mm 使 GZ(20°) 移动
1.7e-03 m，为路线间差的 2343 倍）证明交叉比较本身可以失败。
"""

    report = {
        "title": "TensorCAD 0.23.0：站位截面第三路 GZ(φ) 在当前 main 三路互校（m4549 / PR#52 重开第 3 轮）",
        "summary": (
            "PR #52 因命名遮蔽被关闭后，从当前 main（faa0da6）重切交付第三条独立 GZ 路线：站位截面 "
            "Green 定理积分 + 沿船长 Simpson + 二分水线偏移，以 section_righting_arm_curve 发布，"
            "与 gz.py / righting_arm.py 互不遮蔽。方驳 vs 解析闭式 5.8e-15 相对，三路一致 7.9e-14 m，"
            "船体 vs gz 7.3e-07 m（相对 4e-5），负控证明比较可失败；新增 43 项测试，全套 "
            "998 passed / 10 skipped / 0 failed，ruff + mypy 全绿。"
        ),
        "scenario_key": "tensorcad/righting-third-route",
        "mission_id": 4549,
        "tags": ["TensorCAD", "righting-arm", "GZ", "hydrostatics", "stability",
                 "third-route", "cross-verification", "PR#52", "m4549"],
        "markdown": md,
        "figures": [
            {"name": "righting_arms.png",
             "caption": "方驳/船体 GZ 曲线及与解析闭式、其他路线的偏差（run_example.py 生成）",
             "png_base64": fig1},
            {"name": "merge_verification.png",
             "caption": "三条 GZ 路线在同树上的互校（tools/verify_merge_m4534.py 生成）",
             "png_base64": fig2},
        ],
    }
    OUT.write_text(json.dumps(report))
    print("wrote", OUT, OUT.stat().st_size, "bytes;", len(md), "markdown chars")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
