"""Render the 3-point (10^6->10^7) ladder JSON as markdown tables.

Zero hand-copied numbers: this script is the single source for the tables that
go into ``docs/perf_scale_10m.md`` and into the platform report.  It prints
four blocks and nothing else:

1. per-rung measurements (sub-task, bodies, triangles, seconds, throughput,
   index build, index MB, mesh MB, live RSS, peak RSS),
2. the scaling (complexity) gates with their verdict,
3. the volume accuracy check against the analytic reference,
4. the memory-vs-size summary (mesh footprint, index, live RSS).

Usage::

    python scripts/perf_3pt_report.py benchmarks/perf_scale/ladder_fleet_10m_3pt.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def _fmt(value, digits: int = 3, dash: str = "–") -> str:
    if value is None:
        return dash
    if isinstance(value, float):
        return f"{value:,.{digits}f}"
    return f"{value:,}"


def results_table(data: dict) -> str:
    head = ("| 子任务 | 体数 | 三角面 | 秒 | 吞吐 | 建索引 s | 索引 MB | 网格 MB | 活 RSS MB | 峰值 RSS MB | 状态 |\n"
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n")
    rows = []
    for r in data.get("results", []):
        rows.append(
            "| `{sub}` | {copies} | {tri} | {sec} | {tp} | {build} | {idx} | {mesh} | {rss} | {peak} | {status} |".format(
                sub=r["subtask"], copies=r.get("copies", 1), tri=_fmt(r.get("triangles")),
                sec=_fmt(r.get("seconds"), 3), tp=_fmt(r.get("throughput"), 0),
                build=_fmt(r.get("build_seconds"), 2),
                idx=_fmt((r.get("index_bytes") or 0) / (1024.0 * 1024.0), 1) if r.get("index_bytes") else "–",
                mesh=_fmt(r.get("mesh_mb"), 0), rss=_fmt(r.get("rss_mb"), 0),
                peak=_fmt(r.get("peak_rss_mb"), 0), status=r.get("status"),
            )
        )
    return head + "\n".join(rows)


def gates_table(data: dict) -> str:
    head = ("| 门 | 判定 | 细节 |\n|---|---|---|\n")
    rows = []
    for g in data.get("gates", []):
        mark = "not_run（未检查）" if g["status"] == "not_run" else ("pass" if g["passed"] else "**RED**")
        rows.append(f"| `{g['name']}` | {mark} | {g['detail']} |")
    footer = (
        f"\n\nverdict=**{data.get('verdict')}**，gates_checked={data.get('gates_checked')}，"
        f"gates_not_run={data.get('gates_not_run')}"
    )
    return head + "\n".join(rows) + footer


def volume_table(data: dict) -> str:
    head = "| 体数 | 三角面 | 相对误差 | <3%? |\n|---:|---:|---:|---|\n"
    rows = []
    for c in data.get("volume_checks", []):
        rel = c.get("rel_err")
        rows.append(
            f"| {c.get('copies', 1)} | {_fmt(c.get('triangles'))} | "
            f"{rel * 100:+.4f}% | {'是' if c.get('ok') else '**否**'} |"
        )
    return head + "\n".join(rows) if rows else "（无体积检查行）"


def memory_table(data: dict) -> str:
    head = "| 三角面 | 网格估算 MB | 索引 MB | 活 RSS MB | 峰值 RSS MB |\n|---:|---:|---:|---:|---:|\n"
    rows, seen = [], set()
    for r in sorted(data.get("results", []), key=lambda x: x.get("triangles") or 0):
        if r.get("status") != "measured":
            continue
        key = (r["subtask"], r.get("copies", 1), r.get("triangles"))
        if key in seen:
            continue
        seen.add(key)
        rows.append(
            f"| {_fmt(r.get('triangles'))} | {_fmt(r.get('mesh_mb'), 1)} | "
            f"{_fmt((r.get('index_bytes') or 0) / (1024.0 * 1024.0), 1) if r.get('index_bytes') else '–'} | "
            f"{_fmt(r.get('rss_mb'), 0)} | {_fmt(r.get('peak_rss_mb'), 0)} |"
        )
    return head + "\n".join(rows)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Markdown tables for the 3-point perf ladder")
    ap.add_argument("json", nargs="?", default="benchmarks/perf_scale/ladder_fleet_10m_3pt.json")
    args = ap.parse_args(argv)
    src = Path(args.json)
    if not src.exists():
        raise SystemExit(f"missing {src}")
    data = json.loads(src.read_text(encoding="utf-8"))

    print("## A. 逐档测量\n")
    print(results_table(data))
    print("\n## B. 复杂度/回归门\n")
    print(gates_table(data))
    print("\n## C. 体积精度（vs 解析基准）\n")
    print(volume_table(data))
    print("\n## D. 内存随规模\n")
    print(memory_table(data))
    print(f"\n环境: {data.get('environment')}")
    print(f"命令: {data.get('command')}")
    print(f"总墙钟: {data.get('seconds')} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())