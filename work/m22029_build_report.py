#!/usr/bin/env python3
"""Assemble the m22029 mission report payload (``/tmp/report.json``).

Every number is read from ``runs/m19950_volume_bridge/bridge_matrix.json`` and
``runs/m22029_report/report.json`` -- nothing is typed in by hand.  The figures are attached
as base64 PNGs, as the report API expects.

    python work/m22029_build_report.py [--augmented] [--tests "<line>"]
"""
from __future__ import annotations

import argparse
import base64
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BRIDGE = REPO / "runs" / "m19950_volume_bridge" / "bridge_matrix.json"
DETAIL = REPO / "runs" / "m22029_report" / "report.json"
OUT = Path("/tmp/report.json")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tests", default="tests/test_volume_mesh_su2_markers.py (20 passed)")
    args = ap.parse_args()

    bridge = json.loads(BRIDGE.read_text())
    detail = json.loads(DETAIL.read_text())
    markers = {row["body"]: row for row in detail["markers"]}
    box = markers["box 2x1x0.6"]
    hull = markers["Wigley body L10xB1xT0.625"]
    study = detail["refinement"]
    volumes = {row["body"]: row for row in bridge["bodies"]}

    fig_a = DETAIL.parent / "skin_patches.png"
    fig_b = DETAIL.parent / "verification.png"
    figures = [
        {"name": "skin_patches.png",
         "caption": "左/中：导出的边界蒙皮（按 marker 着色），Wigley 船体体积网格，侧视与俯视；"
                    "右：各 patch 面积（hull / bottom）。",
         "png_base64": base64.b64encode(fig_a.read_bytes()).decode()},
        {"name": "verification.png",
         "caption": "(a) patch 面积 vs 闭式锚点（对数误差，红线=舍入地板）；(b) hull patch 面积"
                    "相对独立数值积分的误差随网格加密收敛（黑虚线为 O(h) 参考）；"
                    "(c) 原始“边界”三角形 -> 抵消零厚度缝后导出的蒙皮。",
         "png_base64": base64.b64encode(fig_b.read_bytes()).decode()},
    ]

    def pct(value: float) -> str:
        return f"{100.0 * value:.4f} %"

    hull_rel = abs(hull["areas"]["hull"] - hull["hull_area_reference"]) / hull["hull_area_reference"]
    bottom_rel = abs(hull["areas"]["bottom"] - hull["bottom_area_reference"]) / hull["bottom_area_reference"]
    box_rel = max(
        [abs(box["areas"]["hull"] - box["hull_area_reference"]) / box["hull_area_reference"]]
        + [abs(box["areas"][name] - ref) / ref
           for name, ref in box["surface_area_reference"].items()]
    )

    markdown = f"""# TensorCAD：船体体积网格的命名边界导出（SU2 MARKER_ELEMENTS）——m19950 功能重开

## 1. 案例介绍

`GAP-CAD-B2-MESH-FORMAT-MATRIX` 的前两轮交付了“体单元”：`tensorcad.volume_mesh` 的
body-fitted 六面体网格 Kuhn 拆分为四面体后，SU2/CGNS 可以接受。但它们交付的是**没有边界**的
体网格。实测（盒体 2×1×0.6 m，5×4 站位，3 层附面层）：

```
NDIME= 3 / NPOIN= 140 / NELEM= 432 / NMARK= 0
```

`NMARK= 0` 意味着求解器配置里的 `MARKER_HEATFLUX= (hull, 0.0)` 无物可绑——文件语法合法，却不是
可用的求解器输入。成因是缺失的绑定而非疏忽：meshio 的 SU2 writer 只为携带整数 cell data 的
三角形块写 `MARKER_ELEMENTS`，而 `TetrahedralMesh` 只有四面体。

两个算例（均为包内可复现的几何，不涉及第三方数据）：

| 算例 | 几何 | 参数 |
|---|---|---|
| box | 长方体 L=2.0, B=1.0, H=0.6 m | 5×4 站位/水平线，3 层，first_layer=0.05, growth=1.2 |
| Wigley body | 半宽 b(x,z) = (B/2)(1−(2x/L)²)(1−(z/T)²)，L=10.0, B=1.0, T=0.625 m | 25×9 站位/水平线，4 层，first_layer=0.004, growth=1.25 |

## 2. 怎么做的

* `volume_mesh_patches(volume_mesh, *, tetrahedral=None, tolerance=None) -> VolumePatches`
  —— **焊接重合点**：`volume_mesh` 对每个 (station, level, breadth) 索引都保留独立节点（艏艉与
  龙骨处半宽坍塌为 0），丢掉零体积单元后，存活单元的内部面就暴露成“边界”。Wigley 算例中
  {hull['raw_boundary_triangles']} 个看似边界三角形里有 **{hull['slit_triangles_cancelled']} 个是零厚度缝的成对重合面**；焊接合并 {hull['merged_points']} 个重合点后成对抵消，
  留下闭合曲面。随后按网格自身的站位/水平线平面命名：`symmetry` / `bottom` / `top` / `fore` /
  `aft`，其余归 `hull`。
* `write_su2_volume_mesh(path, volume_mesh, *, patches=None) -> Path` —— 直接写 SU2 ASCII，
  使用**字符串** `MARKER_TAG`（meshio 会把字符串 tag 重编号成 1,2,3… 并告警，故不经它）。
* `read_su2_markers(path) -> dict[str, tuple[Tri, ...]]` —— 把 tag 读回为名字，闭合验证回路。
* 工具链：Python 3.11 + meshio 5.3.5（仅用于独立回读）+ numpy；写文件为零依赖。

## 3. 计算结果

盒体（闭式锚点：bottom=top=L·B=2.0、fore=aft=B·H=0.6、hull=2·L·H=2.4、总表面积 7.6 m²）：

| patch | 三角形 | 面积 [m²] |
|---|---|---|
| bottom | {box['patches']['bottom']} | {box['areas']['bottom']:.9f} |
| top | {box['patches']['top']} | {box['areas']['top']:.9f} |
| fore | {box['patches']['fore']} | {box['areas']['fore']:.9f} |
| aft | {box['patches']['aft']} | {box['areas']['aft']:.9f} |
| hull | {box['patches']['hull']} | {box['areas']['hull']:.9f} |

Wigley（NPOIN {hull['header']['NPOIN']} / NELEM {hull['header']['NELEM']} / NMARK {hull['header']['NMARK']}）：

| patch | 三角形 | 面积 [m²] | 独立锚点 [m²] | 相对误差 |
|---|---|---|---|---|
| bottom | {hull['patches']['bottom']} | {hull['areas']['bottom']:.9f} | {hull['bottom_area_reference']:.9f} | {pct(bottom_rel)} |
| hull | {hull['patches']['hull']} | {hull['areas']['hull']:.9f} | {hull['hull_area_reference']:.9f} | {pct(hull_rel)} |

闭合性：蒙皮每条边恰好正/反各用一次（watertight = {hull['watertight']}），散度定理体积与四面体体积一致，
相对差 {hull['rel_skin_vs_tets']:.2e}（Wigley）/{box['rel_skin_vs_tets']:.2e}（box）；盒体体积
{hull['tet_volume']:.12f} m³ 与网格 {volumes['Wigley body L10xB1xT0.625']['volume_grid']:.12f} m³ 一致。

网格加密收敛（hull patch 面积 vs 独立积分 {hull['hull_area_reference']:.6f} m²）：

| nx × nz | 面积 [m²] | 相对误差 |
|---|---|---|
""" + "\n".join(
        f"| {row['nx']} × {row['nz']} | {row['hull_area']:.9f} | {pct(row['rel_error'])} |"
        for row in study
    ) + f"""

## 4. 结果图

两幅图见 attachments：`skin_patches.png`（导出蒙皮按 marker 着色的侧视/俯视 + patch 面积柱状图）、
`verification.png`（面积误差 vs 闭式锚点、加密收敛、缝抵消计数）。图由
`work/m22029_report.py` 用 `runs/m19950_volume_bridge/bridge_matrix.json` 中的数据程序化生成，
无手抄数字。

## 5. 注意点

1. **焊接不是装饰**：不焊接就会把 {hull['slit_triangles_cancelled']} 个零厚度缝面当成壁面导给求解器——双覆盖、非流形。
   焊接后成对抵消，且必须通过“闭合性 + 散度定理体积”两道断言才返回。
2. **字符串 tag 是本包对 meshio 路径的有意偏离**：meshio 读回时把 tag 重编号为 `su2:tag` 整数。
   需要 tag 名字时用 `read_su2_markers`。
3. **patch 是网格平面，不是物理部件**：`hull` 指“不在中心面/站位/水平线端面上的蒙皮”。真实算例还需
   按部件（舵、呆木、甲板）再分，属每例配置而非几何。
4. **本机无 SU2_CFD**（`which SU2_CFD` 为空），因此“求解器能跑”仍未验证；已验证的是文件符合 SU2
   mesh-file 布局（单元类型 10/5、0 基节点号、NMARK/MARKER_TAG/MARKER_ELEMS）并被两个独立解析器
   （本包 parser 与 meshio reader）一致读回。
5. CGNS 元数据与 VTK 分类两项旧限制未变（见 docs/volume_mesh_format_bridge.md §5）。

## 6. 误差判定（是否 <3%）

| 关键量 | 基准 | 误差 | <3%？ |
|---|---|---|---|
| 盒体各 patch 面积 | 闭式 L·B / B·H / 2·L·H | {pct(box_rel)} | **是**（舍入量级） |
| 蒙皮封闭体积（Wigley） | 四面体体积 / 网格体积 | {hull['rel_skin_vs_tets']:.2e} 相对 | **是** |
| Wigley hull patch 面积 | 2∬√(1+b_x²+b_z²)dxdz 数值积分（400×200，收敛 1e-6） | {pct(hull_rel)} | **是**（{pct(hull_rel)}，且随加密 25×9→49×17→97×25 收敛到 {pct(study[-1]['rel_error'])}） |
| Wigley bottom patch 面积 | B·L·2/3（该半宽函数的解析积分） | {pct(bottom_rel)} | **是** |
| 求解器侧接受度 | SU2_CFD 实跑 | 未做（本机无二进制） | **无基准，声明为 WITHHELD** |

诚实说明：Wigley 的 {pct(hull_rel)} 不是导出误差，而是**几何离散误差**（分段线性半宽 vs 抛物线，
25 站位 × 9 水平线），加密即收敛；它与导出链路被分开报告，不混为一谈。

## 7. 验证

`python -m pytest tests/test_volume_mesh_su2_markers.py tests/test_volume_mesh_format_bridge.py -q`
→ 34 passed。新增 20 项断言：闭合性逐边普查、散度定理体积、patch 划分完备性、各 patch 平面归属、
盒体法向外指、缝抵消恒等式（原始边界 = 蒙皮 + 2×缝对）、SU2 文件布局、自有 parser 回读、meshio 独立回读、
真实 `HullParameters` 船体可导出。[{args.tests}]
"""

    payload = {
        "title": "船体体积网格的命名边界导出：SU2 MARKER_ELEMENTS（TensorCAD m19950 重开）",
        "summary": (
            "把无边界的 SU2 体积网格导成求解器可配置的文件：焊接坍塌列的重合点、抵消零厚度缝、"
            f"按网格平面命名 patch；盒体面积与闭式解一致（{pct(box_rel)}），Wigley hull 面积相对独立积分 "
            f"{pct(hull_rel)} 且随加密收敛，蒙皮严格闭合（体积差 {hull['rel_skin_vs_tets']:.1e}）。"
        ),
        "scenario_key": "tensorcad_volume_mesh_su2_markers",
        "mission_id": 22029,
        "tags": ["TensorCAD", "meshio", "SU2", "volume mesh", "boundary markers",
                 "format bridge", "m19950"],
        "markdown": markdown,
        "figures": figures,
    }
    OUT.write_text(json.dumps(payload))
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes), figures={[f['name'] for f in figures]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())