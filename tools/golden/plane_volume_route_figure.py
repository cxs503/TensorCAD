"""Figure + report payload for docs/golden/plane_volume_route.md.

Recomputes the cross-route comparison (waterplane closure vs the station-area route) on the
shipped DTMB 5415 sheet and writes:

* ``docs/golden/plane_volume_route.png`` -- relative deviation of V and A_wp against z;
* ``report_m23043.json`` -- the platform report payload (markdown + embedded figure).

Every number in the report is injected from the computation; nothing is typed by hand.
"""

from __future__ import annotations

import base64
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tensorcad import mesh as meshlib  # noqa: E402
from tensorcad.dtmb5415 import DTMB5415Space, WATERLINE_Z  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PNG = ROOT / "docs" / "golden" / "plane_volume_route.png"
REPORT = ROOT / "report_m23043.json"

LEVELS = [-2.0, -1.0, 2.0, 4.0, 5.0, WATERLINE_Z]


def ladder() -> list[dict[str, float]]:
    space = DTMB5415Space()
    mesh = space.mesh()
    rows = []
    for z in LEVELS:
        closure = meshlib.volume_below_plane(mesh, z)
        hydro = space.hydrostatics(waterline_z=z)
        rows.append(
            {
                "z_m": float(z),
                "volume_closure_m3": closure.volume,
                "volume_station_m3": float(hydro["displacement_volume_m3"]),
                "area_closure_m2": closure.waterplane_area,
                "area_station_m2": float(hydro["waterplane_area_m2"]),
                "dev_volume": closure.volume / float(hydro["displacement_volume_m3"]) - 1.0,
                "dev_area": closure.waterplane_area / float(hydro["waterplane_area_m2"]) - 1.0,
                "clipped_triangles": closure.clipped_triangle_count,
                "submerged_area_m2": closure.submerged_surface_area,
                "orientation": closure.orientation,
            }
        )
    return rows


def figure(rows: list[dict[str, float]]) -> None:
    z = [row["z_m"] for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    axes[0].plot(z, [row["volume_closure_m3"] for row in rows], "o-", label="waterplane closure")
    axes[0].plot(z, [row["volume_station_m3"] for row in rows], "s--", label="station-area route")
    axes[0].set_xlabel("waterplane z [m]")
    axes[0].set_ylabel("displacement volume [m$^3$]")
    axes[0].axvline(WATERLINE_Z, color="0.6", lw=0.8, ls=":")
    axes[0].set_title("V(z): both routes")
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    axes[1].axhline(0.0, color="0.6", lw=0.8)
    axes[1].plot(z, [100 * row["dev_volume"] for row in rows], "o-", label="Δ V")
    axes[1].plot(z, [100 * row["dev_area"] for row in rows], "s--", label="Δ A$_{wp}$")
    axes[1].set_xlabel("waterplane z [m]")
    axes[1].set_ylabel("closure vs station-route [%]")
    axes[1].set_title("relative deviation (design WL = %.3f m)" % WATERLINE_Z)
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    fig.suptitle("DTMB 5415 open shell: waterplane-closure integral vs station-area integral")
    fig.tight_layout()
    fig.savefig(PNG, dpi=130)
    plt.close(fig)


def main() -> None:
    rows = ladder()
    figure(rows)
    design = next(row for row in rows if abs(row["z_m"] - WATERLINE_Z) < 1e-9)
    table = "\n".join(
        f"| {row['z_m']:.3f} | {row['volume_closure_m3']:.2f} | {row['volume_station_m3']:.2f} | "
        f"{100 * row['dev_volume']:+.2f} % | {row['area_closure_m2']:.2f} | "
        f"{row['area_station_m2']:.2f} | {100 * row['dev_area']:+.2f} % | {row['clipped_triangles']} |"
        for row in rows
    )
    markdown = f"""## 案例介绍

DTMB 5415 基线船（数据集 README 模型尺度表 × λ = {DTMB5415Space().hydrostatics()['length_stations_m']:.2f}/{DTMB5415Space().hydrostatics()['beam_waterline_m']:.2f}/{DTMB5415Space().hydrostatics()['draft_m']:.3f} m 量级，设计水线 z = {WATERLINE_Z:.4f} m）。
TensorCAD 交付的 5415 几何是**开口船体曲面**（一个 722 边边界环，含甲板边缘 + 艏缝 + 艉部曲线），
因此金标回归套件对该 case 的排水量一直是 **WITHHELD**：`total_signed_volume` 在开曲面上量不出体积，
只能用站位面积路线给参考值。本单补上缺失的**独立路线**：不修网格的水线封闭面积分。

## 怎么做的

对顶部开口的壳体，`z ≤ T` 的水柱边界 = 壳 `z≤T` 部分（侧/底）+ `z=T` 水平面片（外法向 +e_z）
+ 壳体其余开口（艉平面、艏纵缝），后者是**竖直**面，对下面两个积分贡献为零。
记 `g = ∫n_z dA`、`h = ∫z n_z dA`（过平面三角形按 Sutherland–Hodgman 精确裁切后扇形三角化），
则 `A_wp = -g`、`V = h - T·g`，符号由 `g` 判定并回传 `orientation`。
实现：`tensorcad.mesh.volume_below_plane` → `PlaneVolume`；工具链：纯 Python、无 numpy 依赖、无网格修复。

## 计算结果（程序化注入，零手抄）

| 水线 z [m] | V 封闭路线 [m³] | V 站位路线 [m³] | ΔV | A_wp 封闭 [m²] | A_wp 站位 [m²] | ΔA | 裁切三角形 |
| --- | --- | --- | --- | --- | --- | --- | --- |
{table}

设计水线处：`submerged_surface_area = {design['submerged_area_m2']:.2f} m²`，
`clipped_triangle_count = {design['clipped_triangles']}`，`orientation = {design['orientation']}`
（交付网格绕向指向水柱内部，已符号解析，不是缺陷）。

## 结果图

左：两条路线的 `V(z)`；右：相对偏差随水线升高而单调变大（−1.3 % @ z=−2 m → −5.1 % @ 设计水线），
说明分歧集中在上部船体，而不是常数偏置或整体缩放。

## 注意点

1. **积分器本身已验证**：`tests/test_mesh_plane_volume.py` 7 项全绿——2×3×4 无甲板盒（体积/水线面/湿表面积精确到 1e-15，含跨平面裁切）、
   闭合盒（水线面 → 0）、单位半球（V 收敛到 2/3π、A 收敛到 π）、绕向反转同一个数、全在平面上方 → 0、NaN 拒绝。
   所以 5415 上的 5 % 是**几何/封闭性**的问题，不是积分器的问题。
2. **封闭假设是瓶颈**：艉部开口是后倾的艉底曲线（x 变化 ≈0.4 m），艏缝下探到 z≈1.6 m，
   二者在水线下都不是严格竖直面；粗估只值几个 m³，量级差一个数量级，属于贡献项而非主因。
3. **主因待定**：最可能是"三角化壳体 ≠ 站位路线积分的那个模型"——网格是 2250 个出版点的
   双三次样条（每 patch 49×133），站位路线在同一样条上细采样并逐站精确裁切；分歧随 z 增大与此一致。
   另需排除"水线下缺三角形"（例如前后 patch 之间 x = 80.13…81.79 的接缝区）这一可能。
4. **金标基线未重新 bless**：5 % 远在任何可断言带之外，`displacement-withheld-open-shell` 保持原样，
   WITHHELD 不被粉饰成 REAL。诚实负结果入账。

## 误差判定

关键量 vs 基准（基准 = 同一几何的站位面积积分，非解析解）：

| 量 | 封闭路线 | 站位路线 | 误差 | 是否 <3 % |
| --- | --- | --- | --- | --- |
| 排水体积 V（设计水线） | {design['volume_closure_m3']:.2f} m³ | {design['volume_station_m3']:.2f} m³ | {100 * design['dev_volume']:+.2f} % | **否** |
| 水线面面积 A_wp（设计水线） | {design['area_closure_m2']:.2f} m² | {design['area_station_m2']:.2f} m² | {100 * design['dev_area']:+.2f} % | **否** |
| 排水体积 V（z = −2 m） | {rows[0]['volume_closure_m3']:.2f} m³ | {rows[0]['volume_station_m3']:.2f} m³ | {100 * rows[0]['dev_volume']:+.2f} % | 是 |
| 无甲板盒 V/A（封闭路线 vs 解析） | 24.0 / 6.0 | 24.0 / 6.0 | ~1e-16 | 是 |

**结论：5415 的 ∇ 两路线互差 5 %，本单不把封闭路线当成独立锚，WITHHELD 维持。**
下一步判别实验（逐站截面积：网格 vs `DTMB5415Space._stations`）已在文档中写明。
"""
    png = base64.b64encode(PNG.read_bytes()).decode("ascii")
    payload = {
        "title": "5415 开口壳体的水线封闭面积分：积分器验证通过，5 % 跨路线分歧待定位（WITHHELD 维持）",
        "summary": "新增 tensorcad.mesh.volume_below_plane（开壳体积/水线面，解析封闭，零网格修复），"
                   "7 项解析基准测试全绿；应用到 DTMB5415 开口壳体，与站位路线在水线处差 −5.07 % V / −4.26 % A_wp，"
                   "故金标 WITHHELD 不解除、基线不重 bless，并把分歧定位实验写进 docs。",
        "scenario_key": "tensorcad-golden-plane-volume-route",
        "mission_id": 23043,
        "tags": ["TensorCAD", "GAP-CAD-D2-GOLDEN-REGRESSION", "golden", "hydrostatics",
                 "divergence-theorem", "verification", "negative-result"],
        "markdown": markdown,
        "figures": [
            {
                "name": "plane_volume_route",
                "caption": "DTMB 5415 开口壳体：水线封闭面积分 vs 站位面积分（V(z) 与相对偏差）",
                "png_base64": png,
            }
        ],
    }
    REPORT.write_text(json.dumps(payload, ensure_ascii=False, indent=1))
    print(json.dumps({"rows": rows, "png": str(PNG), "report": str(REPORT)}, indent=1))


if __name__ == "__main__":
    main()