"""Figures + progress report for mission m21617 (GAP-CAD-D4-PERF-SCALE evidence bundle)."""
from __future__ import annotations

import base64
import io
import json
import os
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/root/tp-agent/work/m21617")
sys.path.insert(0, str(ROOT))
from tensorcad.perf_evidence import (  # noqa: E402
    load_artifact, render_budget_table, render_mesh_table, render_voxel_table,
)

MESH = load_artifact(ROOT / "benchmarks/perf_scale/artifacts/ladder_mesh.json")
VOXEL = load_artifact(ROOT / "benchmarks/perf_scale/artifacts/ladder_voxel.json")
FIGS = ROOT / "benchmarks/perf_scale/figures"
FIGS.mkdir(parents=True, exist_ok=True)


def series(sub: str):
    pts = [(r["triangles"], r["seconds"]) for r in MESH["results"]
           if r["subtask"] == sub and r["status"] == "measured" and r["seconds"]]
    return sorted(pts)


# ---- figure 1: wall clock vs triangles (log-log) ----
fig, ax = plt.subplots(figsize=(7.2, 4.6))
styles = {"load": ("o", "#1f77b4"), "point_locate": ("s", "#d62728"),
          "point_locate_accel": ("^", "#2ca02c"), "volume": ("d", "#9467bd")}
for name, (mk, col) in styles.items():
    pts = series(name)
    if not pts:
        continue
    xs, ys = zip(*pts)
    ax.plot(xs, ys, marker=mk, color=col, label=name)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("triangles"); ax.set_ylabel("wall clock [s]")
ax.set_title("TensorCAD 10^6-band ladder (mission m21617, measured)")
ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout()
p1 = FIGS / "scaling_m21617.png"; fig.savefig(p1, dpi=130); plt.close(fig)

# ---- figure 2: voxelisation ----
fig, ax = plt.subplots(figsize=(6.0, 4.2))
vr = sorted((r for r in VOXEL["results"] if r["status"] == "measured"),
            key=lambda r: r["grid_cells"])
xs = [r["grid_cells"] for r in vr]; ys = [r["seconds"] for r in vr]
ax.plot(xs, ys, marker="o", color="#ff7f0e", label="voxelize (measured)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("voxel cells"); ax.set_ylabel("wall clock [s]")
ax.set_title("voxelisation ladder (32^3 .. 128^3)")
ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout()
p2 = FIGS / "voxelize_m21617.png"; fig.savefig(p2, dpi=130); plt.close(fig)

# ---- figure 3: error/invariant panel (doc-vs-artifact + volume agreement) ----
fig, ax = plt.subplots(figsize=(6.4, 3.4))
ax.axis("off")
ax.text(0.01, 0.95, "Checked invariants (this increment)", fontsize=11, weight="bold", va="top")
ax.text(0.01, 0.80,
        "doc tables vs committed artifacts: 0 drift (3 marked tables, exact match)\n"
        "cross-machine baseline: refused (status=not_run, never pass)\n"
        "accel vs brute-force point location: identical hit counts, same grid\n"
        "volume vs analytic 7.512698 m3 (previous run): -0.0030% @4.0e6 tri",
        fontsize=8.5, va="top", family="monospace")
fig.tight_layout()
p3 = FIGS / "invariants_m21617.png"; fig.savefig(p3, dpi=130); plt.close(fig)


def b64(path: Path) -> str:
    return base64.b64encode(path.read_bytes()).decode()


# ---- headline numbers, programmatic ----
accel = dict((t, s) for t, s in series("point_locate_accel"))
load = dict((t, s) for t, s in series("load"))
big = max(load)
mesh_rows = len(MESH["results"])
vox_rows = len(VOXEL["results"])
env = MESH["environment"]

markdown = f"""## 1. 案例介绍

TensorCAD（零依赖几何内核）的大几何性能基准锚点 `GAP-CAD-D4-PERF-SCALE`：参数化船体
（L=10 m, B=2 m, T=0.8 m）网格阶梯，甲板到 10⁶ 三角级；子任务 = 加载 / 点定位（暴力 vs 均匀网格+DDA 加速）/
体素化 / 体积积分。本单不是再测一遍，而是**把已合并的测量结论变成可复核的证据包**。

实测环境（来自工件 `environment` 字段，非手抄）：`{env.get('machine')}` / `{env.get('platform')}` /
cpu_count={env.get('cpu_count')} / python {env.get('python')} / TensorCAD {env.get('tensorcad')}。

## 2. 怎么做的

1. **审计**（对抗式复核）：文档声称"全部来自 `artifacts/perf_scale/full_ladder.json`，程序化读取"，
   但该目录被合并门禁滤除、仓库内不存在 → §3 的 44 条测量在干净检出上不可审计；表格数字是手抄进
   markdown 的，无任何机制阻止与工件漂移；`evaluate_gates()` 记录基线 `environment` 却从不核对，
   换机器后会拿别人的秒数判 RED/判绿。
2. **工件入库**：`benchmarks/perf_scale/artifacts/{{ladder_mesh,ladder_voxel}}.json` 受版本控制
   （{mesh_rows} 条网格测量 + {vox_rows} 条体素测量，最大档 {big:,} 三角面）。
3. **表格程序化注入**：`tensorcad/perf_evidence.py` 渲染三张表（网格 / 体素 / 60 s 预算等价规模），
   文档里夹在 `<!-- perf-evidence:... -->` 标记之间；`scripts/perf_scale_doc.py --write/--check`。
4. **基线环境闸门**：`environment_compat()`；`machine/platform` 不一致 → `regression_environment`
   门 `status=not_run`（不是 pass），不产生任何 `regression_*` 绿门。
5. 验证：`pytest tests/test_perf_scale_evidence.py tests/test_perf_benchmark.py` → **57 passed**。

## 3. 计算结果（全部由工件渲染，零手抄）

### 3.1 网格子任务

{render_mesh_table(MESH)}

### 3.2 体素化

{render_voxel_table(VOXEL)}

### 3.3 60 s 预算等价规模（规则：实测速率 × 60 s）

{render_budget_table(MESH, VOXEL)}

资产峰值：最大档 {big:,} 三角面，`load` {load.get(big, float('nan')):.1f} s；
`point_locate_accel` 同档 {accel.get(big, float('nan')):.3f} s。

## 4. 结果图

见 figures：scaling_m21617.png（墙钟 vs 三角面 log-log，四条通道）、voxelize_m21617.png
（体素阶梯）、invariants_m21617.png（本单校验的不变量清单）。

## 5. 注意点 / 局限

- **本机有界阶梯**：`21×9 … 1281×321`（最高 {big:,} 三角面）在本机实测；10⁷ 带本次**未重测**
  （§8.3 的 `cj3081` 多体实测仍是唯一 10⁷ 数据）。
- 递交的全阶梯作业 `cj3563`（网格）/`cj3564`（体素）在本单窗口内始终 `queued`（集群约 35 个作业在前），
  未取得回执 → 未覆盖原 n315 全阶梯（4.0×10⁶ / 384³）。到件后 `perf_scale_doc.py --write` 即可换表。
- 跨机不可逐档比：本机 load1≈70/128，n315 当时≈40。
- `peak_rss_mb` 是进程高水位（`ru_maxrss`），不可逐档读；逐档看 `rss_mb`/`mesh_mb`。

## 6. 误差判定（关键量 vs 基准，逐项）

| 量 | 基准 | 实测 | 误差 | <3%? |
|---|---|---|---|---|
| 文档表格 vs 已提交工件（3 张表，逐字符） | `perf_evidence.render_block()` | 完全一致 | 0 | 是（精确相等） |
| 手改一个数字的可检出性 | 注入漂移 | `--check` 报 DRIFT、测试 RED | 检出 | 是（反证测试） |
| 跨架构基线比较 | 不得产生 RED/pass | `regression_environment` `not_run` | 拒绝比较 | 是（反证测试） |
| 体积 vs 解析 7.512698 m³ @4.0×10⁶ tri（上一轮实测） | 解析值 | 7.512474 | −0.0030% | 是 |
| `point_locate_accel` 命中数 vs 暴力参考 | 同网格暴力解 | 逐档完全一致 | 0 | 是 |

**结论**：本增量不含新的物理基准量；它校验的是**证据约束**（文档 == 工件、跨机基线拒绝比较），
这三项都是精确判据、误差为 0。物理精度仍沿用上一轮实测（体积 −0.0030%，远小于 3%）。
诚实边界：10⁷ 带未重测、集群回执未到，见 §5。
"""

report = {
    "title": "TensorCAD 大几何性能基准：证据包可复核化（表格==工件、跨机基线拒绝比较）",
    "summary": (
        "GAP-CAD-D4-PERF-SCALE 增量：把 10^6 带性能表从手抄改为工件渲染（提交 artifacts，"
        "docs 表格由 perf_evidence 注入并由测试守漂移），并为 2× 回归门加基线环境闸门"
        "（跨架构基线→not_run 而非假绿）。测得 1,640,958 三角面阶梯；57 项测试通过。"
    ),
    "scenario_key": "cad.perf_scale.evidence",
    "mission_id": 21617,
    "tags": ["TensorCAD", "performance", "benchmark", "reproducibility", "GAP-CAD-D4-PERF-SCALE"],
    "markdown": markdown,
    "figures": [
        {"name": "scaling_m21617.png",
         "caption": "墙钟 vs 三角面（log-log）：load / 暴力点定位 / 网格+DDA 点定位 / 体积积分，本单实测阶梯至 1.64e6 三角面",
         "png_base64": b64(p1)},
        {"name": "voxelize_m21617.png",
         "caption": "体素化阶梯 32³–128³ 的墙钟随体素数的缩放（实测）",
         "png_base64": b64(p2)},
        {"name": "invariants_m21617.png",
         "caption": "本单校验的不变量：文档-工件零漂移、跨机基线拒绝比较、加速点定位与暴力解一致",
         "png_base64": b64(p3)},
    ],
}

base = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1").rstrip("/")
user = os.environ["TENSOR_AGENT_USER"]
pw = os.environ["TENSOR_AGENT_PASS"]
login = subprocess.run(["curl", "-s", "-X", "POST", f"{base}/auth/login",
                        "-H", "Content-Type: application/json",
                        "-d", json.dumps({"username": user, "password": pw})],
                       capture_output=True, text=True)
tok = json.loads(login.stdout).get("access_token") or json.loads(login.stdout).get("token")
if not tok:
    print("LOGIN FAILED", login.stdout[:400]); raise SystemExit(1)
Path("/tmp/report.json").write_text(json.dumps(report))
post = subprocess.run(["curl", "-s", "-X", "POST", f"{base}/tensor/reports",
                       "-H", f"Authorization: Bearer {tok}",
                       "-H", "Content-Type: application/json",
                       "-d", "@/tmp/report.json"], capture_output=True, text=True)
print("POST", post.stdout[:400])
print("max tri", big, "load s", load.get(big), "accel s", accel.get(big))