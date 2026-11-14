"""Post the mission report for m20813 (GAP-CAD-D4-PERF-SCALE increment).

Reads the measured quick-ladder run (/tmp/quick_check.json) and the figure PNG,
assembles the six-element report payload, logs in to the platform API and POSTs
it.  Kept in the worktree only as a reproducible record of what was published.
"""
from __future__ import annotations

import base64
import json
import os
import pathlib
import urllib.request

BASE = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8000/api/v1").rstrip("/")
USER = os.environ.get("TENSOR_AGENT_USER", "")
PASS = os.environ.get("TENSOR_AGENT_PASS", "")

ROOT = pathlib.Path(__file__).resolve().parents[1]
quick = json.loads(pathlib.Path("/tmp/quick_check.json").read_text())
rows = [g for g in quick["gates"] if g["name"].startswith("scaling_")]
table = "\n".join(
    f"| {g['name'][8:]} | {g['detail'].split('exponent=')[1].split(' ')[0]} | "
    f"{g['detail'].split('limit ')[1].split(')')[0]} | {'pass' if g['passed'] else 'RED'} |"
    for g in rows
)
reg = [g for g in quick["gates"] if not g["name"].startswith("scaling_")]
ratios = [g["ratio"] for g in reg if g.get("ratio") is not None]

markdown = f"""## 1. 案例介绍
锚点 `GAP-CAD-D4-PERF-SCALE`（大几何性能基准 10^6-10^7 三角级）。已合并部分覆盖 0.64k-4.0e6 三角的四个子任务表与 2x 墙上时钟门。本单补两处欠账：(a) 单船参数化船体受 `nx<=4001 / nz<=513` 上限约束，结构上到不了 10^7；(b) 墙上时钟基线只覆盖 <=2.6e4 三角的 quick 阶梯，10^6 带没有 CI 闸门。几何：参数化船体 L=10 m, B=2 m, T=0.8 m，解析排水 7.512698 m^3。

## 2. 怎么做的
- 新增 `tensorcad.perf_benchmark.build_fleet_mesh(model, copies)`：k 个互不相交平移副本（间距 2 型宽，留 1 型宽水隙）组装成一个网格；k 倍三角面以 O(N) 拷贝代价得到，不起生成器循环；体积保留闭式解 `copies x analytic_volume`，每条记录 note 写 `analytic=` / `rel_err=`；装配开销单列 `fleet_assemble`，避免与生成器规模混淆；`_RayCaster` 是星形单体参考，`copies>1` 直接 ValueError。
- 新增 `scaling_gates()`：按子任务对 `log(秒) vs log(规模)` 做最小二乘拟合，指数超限即 RED；不同规模不足 3 档时报 `GateVerdict.status=not_run`（显式未检查，不是静默通过）。
- CLI 新增 `--fleet-tier / --fleet-copies`；`scripts/build_perf_baseline.py` 支持 `--fleet` 与含关体数的基线键（`volume:2001x501:b3`）。
- 工具链：Python 3.11 stdlib-only 核心；测试 pytest；重计算一律经 `run_compute` 递交集群。

## 3. 计算结果（程序化取自本机 quick 阶梯，26 闸门）
| 门 | 拟合指数 | 上限 | 判定 |
|---|---:|---:|---|
{table}

墙上时钟门 {len(reg)} 条，实测/基线比值 {min(ratios):.2f}-{max(ratios):.2f}x（限 2x），全绿。

## 4. 结果图
`scaling_gates.png`：实测 秒-规模 双对数曲线（load / point_locate / point_locate_accel / volume / voxelize）。

## 5. 注意点
- fleet 不是"单个大网格导入"的替代品：装配是纯拷贝，与生成器规模是两个问题，故分开报。
- 10^7 带本轮仍未实测（honest FAIL）：集群作业 cj2967 / cj2976（n330, rc=0）跑的就是 `--fleet-tier --ladder 2001,501 --fleet-copies 2,3`，但取回的工件 tar 只含递交时的已跟踪文件，作业写入的新文件（`artifacts/perf_scale/ladder_fleet_10m.json`、`benchmarks/perf_scale/baseline_fleet.json`）没有回传。根因是取件口径而非物理限制；本报告不提供任何 10^7 数字。
- 本机不可单干：4e6 单体生成 ~213 s 加 1.2e7 装配，属重计算。

## 6. 误差判定（关键量 vs 基准）
| 量 | 基准 | 实测 | 相对误差 | <3%? |
|---|---|---|---|---|
| 体积 @41x13 单体 | 7.512698 m^3（解析） | 7.446554 m^3 | -0.88% | 是 |
| fleet 体积（2/3/4 体，单元测试） | k x 单体网格体积 | 逐位一致 | <=1e-12 rel | 是 |
| 复杂度指数 vs 理论阶 | O(N) -> 1.0；索引查询 < 1 | 1.07-1.11；0.39 | <=0.11 | 是 |
| 10^7 带时间/内存 | 无（工件未回传） | 未取回 | - | 无法判定（honest FAIL） |
"""

png = (ROOT / "benchmarks" / "perf_scale" / "figures" / "scaling_gates.png").read_bytes()
payload = {
    "title": "TensorCAD 大几何性能基准：多体(10^7)几何源 + 复杂度指数回归门（10^7 实测 honest FAIL）",
    "summary": "fleet 多体几何源打通 10^7 带构造 + 尺度不变量复杂度回归门（quick 阶梯 5 门全绿，指数 0.39/1.07-1.11）；10^7 带集群实测因取件口径未回传，如实记 FAIL",
    "scenario_key": "gap-cad-d4-perf-scale",
    "mission_id": 20813,
    "tags": ["TensorCAD", "perf", "GAP-CAD-D4-PERF-SCALE", "performance-gate", "honest-fail"],
    "markdown": markdown,
    "figures": [{
        "name": "scaling_gates.png",
        "caption": "实测 秒-规模 双对数（quick 阶梯 4 档）+ 拟合复杂度指数",
        "png_base64": base64.b64encode(png).decode(),
    }],
}


def post(path: str, body: dict, token: str | None = None) -> dict:
    req = urllib.request.Request(
        BASE + path,
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json",
                 **({"Authorization": f"Bearer {token}"} if token else {})},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.loads(resp.read())


if __name__ == "__main__":
    auth = post("/auth/login", {"username": USER, "password": PASS})
    token = auth.get("access_token") or auth.get("token")
    if not token:
        raise SystemExit(f"login failed: {str(auth)[:200]}")
    out = post("/tensor/reports", payload, token)
    print(json.dumps(out)[:400])