"""Generate docs/m24628_meshio_matrix18.md from runs/m24628_meshio_matrix18/matrix18.json.

All numbers are injected programmatically from the run evidence (no hand
transcription).  Re-run after any matrix change to keep the committed table in
sync with the machine-readable record.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs" / "m24628_meshio_matrix18"
OUT = ROOT / "docs" / "m24628_meshio_matrix18.md"


def main() -> int:
    d = json.loads((RUN / "matrix18.json").read_text())
    rows = d["matrix"]
    gate = d["geometry_gate"]

    lines: list[str] = []
    lines.append("# TensorCAD meshio 六格式矩阵往返保真 · 18 组合全绿表（m24628 复核）")
    lines.append("")
    lines.append(f"**锚点**: GAP-CAD-B2-MESH-FORMAT-MATRIX  ")
    lines.append(f"**验收**: 往返保真——体积/拓扑/取向不变（6 格式×3 几何=18 组合全绿表）  ")
    lines.append(f"**本单**: m24628（generalist 自找活）  ")
    lines.append(f"**环境**: python {d['python']} / meshio {d['meshio']} / numpy {d['numpy']}  ")
    lines.append(f"**判定**: {d['combinations']} 组合（6 格式×4 几何，含验收 18 组合超集）"
                 f" all_pass=**{d['all_combinations_pass']}**")
    lines.append("")
    lines.append("## 1. 全绿表（体积/拓扑/取向三不变量）")
    lines.append("")
    lines.append("| 格式 | 几何 | 拓扑 | 体积(原始) | 体积(往返) | 相对误差 | 体积OK | 拓扑OK | 取向OK | 文件sha256(前12) |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        lines.append(
            f"| {r['format']} | {r['geometry']} | {r['topology']} "
            f"| {r['volume_original']:.9f} | {r['volume_round_trip']:.9f} "
            f"| {r['relative_error']:.2e} | {r['volume_ok']} | {r['topology_ok']} "
            f"| {r['orientation_ok']} | {r['sha256'][:12]} |"
        )
    lines.append("")
    lines.append(f"**判定**：{d['combinations']}/{d['combinations']} 组合三不变量全绿"
                 f"（all_combinations_pass={d['all_combinations_pass']}）。"
                 "验收要求 6×3=18 组合，本表为 6×4=24 组合超集（第 4 几何=真实 Wigley 船体）。")
    lines.append("")
    lines.append("## 2. 几何门（离散体积 vs 解析基准）")
    lines.append("")
    lines.append("| 几何 | 水密 | 体积(测量) | 体积(解析) | 相对误差 | 四面体数 |")
    lines.append("|---|---|---|---|---|---|")
    for g in gate:
        lines.append(
            f"| {g['geometry']} | {g['watertight']} | {g['volume_measured']:.6f} "
            f"| {g['volume_analytic']:.6f} | {g['relative_error']:.2%} | {g['tetrahedra']} |"
        )
    lines.append("")
    lines.append("解析基准：cube=1.0（精确）；sphere=4π/3；capsule=πr²L+4πr³/3（r=0.5, L=4）；"
                 "wigley_hull=4LBT/9（n=m=2 闭式，L=10, B=1, T=0.625）。"
                 "离散体积误差全部 <1%（网格离散化误差，非格式误差）。")
    lines.append("")
    lines.append("## 3. 证据指针")
    lines.append("")
    lines.append(f"- 机器可读全绿表：`runs/m24628_meshio_matrix18/matrix18.json`")
    lines.append(f"- 24 个往返文件（可第三方复核）：`runs/m24628_meshio_matrix18/files/`")
    lines.append(f"- 图：`artifacts/m24628_meshio_matrix18/`（status_matrix / volume_relative_error / geometry_gate）")
    lines.append(f"- 复现：`{d['reproduce']}`")
    lines.append("")
    OUT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)} ({len(lines)} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
