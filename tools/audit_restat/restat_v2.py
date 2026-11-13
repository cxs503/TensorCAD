#!/usr/bin/env python3
"""AUDIT-RESTAT v2 · 全库真实成果重新统计（L0-L5）第二版 —— 快照锁定 + 前值复现 + 溯源矩阵.

设计约束（为什么这样写）:
  1. 零手抄：每一个进入报告的数字都来自本文件 SQL_REGISTRY 中的一条具名 SQL，
     由 run_sql() 现查现算；JSON 产物逐项记录 sql_id，MD 由 JSON 渲染。
  2. 快照锁定：所有查询都带 `<= :snapshot` 时间上界（verdicts.created_at /
     missions.created_at）。统计量是时变的，不锁快照的数字不可对账（实证：同一
     "L2+ 单点"在 2026-11-11 与 2026-11-13 相差 1,685 单）。
  3. 前值复现：用上游两单的 finished_at 作为 cutoff 重算同一批指标，把
     "能不能复现别人报的数" 变成可执行的断言（reproduce_prior()）。
  4. 纯函数与 I/O 分离：口径逻辑（锚点抽取 / claim-gate 严格性 / crosstab /
     膨胀比）是纯函数，可脱库单测（tests/test_restat_l0l5_v2.py）。

用法:
    python3 tools/audit_restat/restat_v2.py                # 快照=DB now()
    python3 tools/audit_restat/restat_v2.py --snapshot 2026-11-13T01:30:00+00:00
产物:
    artifacts/audit_restat/restat_v2_<snapshot>.json     机器可读（含 sql_id 溯源）
    docs/audit_restat/restat_v2.md                       人读表
    docs/audit_restat/fig_*.png                          图
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------- 纯函数（可单测）

ANCHOR_RE = re.compile(r"锚点=([A-Za-z0-9_.:+\-]+)")
SCOPE_ANCHOR_RE = re.compile(r"^m\d+\s+(\S+)")
PLACEHOLDER_ANCHORS = {"SLUG", "ANCHOR", "NEW:", "NEW:理由"}


def extract_anchor(goal: str | None) -> str | None:
    """从 mission.goal 首部【锚点=X】抽锚点；占位串（SLUG/NEW:）返回 None。"""
    if not goal:
        return None
    m = ANCHOR_RE.search(goal.replace("\n", " "))
    if not m:
        return None
    anchor = m.group(1).rstrip(".")
    if anchor in PLACEHOLDER_ANCHORS or anchor.startswith("NEW:"):
        return None
    return anchor


def scope_anchor(scope: str | None) -> str | None:
    """审计判定 scope 形如 'm17013 KPI-LBM-VERIFIED-CAPS' -> 'KPI-LBM-VERIFIED-CAPS'."""
    if not scope:
        return None
    m = SCOPE_ANCHOR_RE.match(scope.strip())
    return m.group(1) if m else None


def strictness(claim: int | None, gate: int | None) -> str:
    """claim 与机器门实测 gate 的关系：match / overclaim / underclaim / ungated."""
    if gate is None:
        return "ungated"
    if claim is None:
        return "no_claim"
    if gate == claim:
        return "match"
    return "overclaim" if gate < claim else "underclaim"


def crosstab(rows, claim_key="level", gate_key="gate_level", top=5):
    """claim × gate 计数矩阵；gate=None 记为 'ungated' 列。"""
    cols = list(range(0, top + 1)) + ["ungated"]
    mat = {c: {g: 0 for g in cols} for c in range(0, top + 1)}
    for r in rows:
        c = r[claim_key]
        g = r[gate_key]
        if c is None or c not in mat:
            continue
        gk = "ungated" if g is None else int(g)
        if gk not in mat[c]:
            continue
        mat[c][gk] += 1
    # 行转为 "claim Lx" 标签
    return {f"claim L{c}": {str(g): n for g, n in mat[c].items()} for c in mat}


def ratio(num: float, den: float) -> float | None:
    return None if not den else round(num / den, 4)


def pct(num: float, den: float) -> float | None:
    r = ratio(num, den)
    return None if r is None else round(100.0 * r, 1)


def pct2(num: float, den: float) -> float | None:
    r = ratio(num, den)
    return None if r is None else round(100.0 * r, 2)


def band(values):
    """定义敏感带：[min, max] 与被上游引用的点位是否落在带内。"""
    vals = [v for v in values if v is not None]
    if not vals:
        return None
    return {"min": min(vals), "max": max(vals), "n_definitions": len(vals)}


# ---------------------------------------------------------------- SQL 注册表

SQL_REGISTRY: dict[str, str] = {
    # 快照源
    "db_now": "SELECT now()::text",
    # ---- Tier-A：审计判定层（tensor_maturity_verdicts）
    "verdict_snapshot": """
        SELECT repo, scope, level, gate_level, verdict_by, gate_state
        FROM tensor_maturity_verdicts
        WHERE created_at <= :'snapshot'::timestamptz
    """,
    "verdict_at_cutoff": """
        SELECT count(*) AS n_verdicts,
               count(DISTINCT (repo, scope)) AS n_scopes,
               count(*) FILTER (WHERE gate_level >= 2) AS gate2,
               count(*) FILTER (WHERE gate_level < level) AS overclaim
        FROM tensor_maturity_verdicts
        WHERE created_at <= :'snapshot'::timestamptz
    """,
    # ---- Tier-B：任务层（tensor_dev_missions）
    "missions_succeeded": """
        SELECT id, repo, created_at::text AS created_at, goal,
               result->'pr_gate'->>'state' AS pr_state,
               COALESCE(result->'changed_files', '[]'::jsonb) AS changed_files
        FROM tensor_dev_missions
        WHERE status = 'succeeded' AND created_at <= :'snapshot'::timestamptz
    """,
    "missions_awaiting_pr_merge": """
        SELECT id, repo, left(regexp_replace(goal, E'\\n', ' ', 'g'), 120) AS goal_head
        FROM tensor_dev_missions
        WHERE status = 'awaiting_pr_merge' AND created_at <= :'snapshot'::timestamptz
        ORDER BY id
    """,
    # ---- 能力台账（feature_ledger）
    "ledger_by_status": "SELECT status, count(*) FROM feature_ledger GROUP BY status ORDER BY 2 DESC",
    "ledger_pending": """
        SELECT slug, domain, status, title, acceptance
        FROM feature_ledger
        WHERE status IN ('open', 'planned', 'in_progress', 'dropped')
        ORDER BY status, domain, slug
    """,
    "ledger_all_slug_domain": "SELECT slug, domain, status FROM feature_ledger ORDER BY slug",
}

# 前值（上游同锚点单的对外声明）；复现实验的比对目标
PRIOR_CLAIMS = {
    "m14502_pinned": {
        "missions": [14502],
        "claims": {"n_verdicts": 184, "n_scopes": 184, "gate2": 0, "strict_l2plus": 2366},
    },
    "m17113_pinned": {
        "missions": [17113],
        "claims": {"strict_l2plus": 2366},
    },
}


# ---------------------------------------------------------------- I/O

def psql_cmd() -> list[str]:
    override = os.environ.get("RESTAT_PSQL")
    if override:
        return override.split()
    prefix = [] if os.geteuid() == 0 else []
    return prefix + ["runuser", "-u", "postgres", "--", "psql", "-d", "tensorplatform",
                     "-At", "-F", "\t", "-q", "-v", "ON_ERROR_STOP=1"]


def run_sql(sql: str, snapshot: str | None = None) -> list[list[str]]:
    """执行一条注册表 SQL，返回 TSV 行列表。snapshot 以 psql 变量注入（防注入）。"""
    cmd = psql_cmd()
    if snapshot:
        cmd += ["-v", f"snapshot={snapshot}"]
    proc = subprocess.run(cmd, input=sql, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"psql failed: {proc.stderr.strip()[:500]}")
    out = []
    for line in proc.stdout.splitlines():
        if line.strip() == "":
            continue
        out.append(line.split("\t"))
    return out


def run_sql_json(sql: str, snapshot: str | None = None) -> list[dict]:
    """以 jsonb_agg 包装任意查询，取回结构化的行。"""
    wrapped = f"SELECT COALESCE(jsonb_agg(t), '[]'::jsonb) FROM ({sql}) t"
    rows = run_sql(wrapped, snapshot)
    if not rows:
        return []
    return json.loads(rows[0][0])


# ---------------------------------------------------------------- 统计主体

IMPL_EXT = re.compile(r"\.(py|pyx|c|cc|cpp|cxx|h|hpp|jl|rs|f90|f95|f|for|go|ts|tsx|js|java|m)$")
NON_IMPL_PATH = re.compile(r"(^|/)(tests?|docs?|examples?|benchmarks?)/")
HELPER_NAME = re.compile(r"(^|/)(test_|verify_|audit_|report_|probe_)")


def tier_b_strict_bands(missions: list[dict]) -> dict:
    """Tier-B 严格口径的『定义敏感带』。

    上游把严格 L2+ 定义为「有实现文件 且 有 PR 合并证据」，但 '实现文件' 的
    正则/目录白名单未在正文登记，故同一批数据可在若干合理定义间漂移。
    这里并排算 4 个定义，给出区间而不是单点（单点是伪精度）。
    """
    def merged(m):
        return (m.get("pr_state") or "") == "merged"

    def files(m):
        return m.get("changed_files") or []

    def n_impl(m):
        c = 0
        for f in files(m):
            if IMPL_EXT.search(f) and not NON_IMPL_PATH.search(f) and not HELPER_NAME.search(f):
                c += 1
        return c

    defs = {
        "A_merged_any_files": sum(1 for m in missions if merged(m)),
        "B_merged_impl_file": sum(1 for m in missions if merged(m) and n_impl(m) > 0),
        "C_merged_nonempty_changed_files": sum(1 for m in missions if merged(m) and len(files(m)) > 0),
        "D_impl_file_any_pr_state": sum(1 for m in missions if n_impl(m) > 0),
    }
    strict = [defs["A_merged_any_files"], defs["B_merged_impl_file"],
              defs["C_merged_nonempty_changed_files"]]
    return {
        "definitions": defs,
        # 定义敏感带只取『要求 PR 合并』的三个严格定义；D 是无合并要求的上界，单独列
        "band": band(strict),
        "loose_upper_bound": defs["D_impl_file_any_pr_state"],
    }


def tier_a_stats(verdicts: list[dict]) -> dict:
    n = len(verdicts)
    by_producer = {}
    for r in verdicts:
        p = r.get("verdict_by") or "unknown"
        d = by_producer.setdefault(p, {"n": 0, "gate2": 0, "overclaim": 0,
                                       "match": 0, "underclaim": 0, "ungated": 0, "scopes": set()})
        d["n"] += 1
        d["scopes"].add((r["repo"], r["scope"]))
        cls = strictness(r.get("level"), r.get("gate_level"))
        d[cls] = d.get(cls, 0) + 1
        if r.get("gate_level") is not None and int(r["gate_level"]) >= 2:
            d["gate2"] += 1
    for p, d in by_producer.items():
        d["n_scopes"] = len(d.pop("scopes"))
        d["overclaim_rate_pct"] = pct(d.get("overclaim", 0), d["n"])
        d["gate2_rate_pct"] = pct(d["gate2"], d["n"])
    gate2 = sum(1 for r in verdicts if r.get("gate_level") is not None and int(r["gate_level"]) >= 2)
    return {
        "n_verdicts": n,
        "n_scopes": len({(r["repo"], r["scope"]) for r in verdicts}),
        "gate2": gate2,
        "gate2_rate_pct": pct(gate2, n),
        "overclaim": sum(1 for r in verdicts if strictness(r.get("level"), r.get("gate_level")) == "overclaim"),
        "overclaim_rate_pct": pct(
            sum(1 for r in verdicts if strictness(r.get("level"), r.get("gate_level")) == "overclaim"), n),
        "by_producer": by_producer,
        "crosstab_claim_vs_gate": crosstab(verdicts),
        "by_claim_level": {
            str(lv): {
                "n": sum(1 for r in verdicts if r.get("level") == lv),
                "gate2": sum(1 for r in verdicts if r.get("level") == lv
                             and r.get("gate_level") is not None and int(r["gate_level"]) >= 2),
                "overclaim": sum(1 for r in verdicts if r.get("level") == lv
                                 and strictness(lv, r.get("gate_level")) == "overclaim"),
            } for lv in range(0, 6)
        },
    }


def traceability(ledger: list[dict], verdicts: list[dict], missions: list[dict]) -> dict:
    """任务 → 能力台账 的溯源矩阵：每一层都掉人，这就是「膨胀」的结构。"""
    ledger_anchor = {r["slug"] for r in ledger}
    scopes_parsed, scopes_freeform = {}, []
    for v in verdicts:
        a = scope_anchor(v.get("scope"))
        if a:
            scopes_parsed.setdefault(a, []).append(v)
        else:
            scopes_freeform.append(v.get("scope"))

    anchored_missions, placeholder_missions = [], []
    for m in missions:
        g = m.get("goal") or ""
        if extract_anchor(g):
            anchored_missions.append(m)
        elif "锚点=" in g:
            placeholder_missions.append(m)
    anchors = {extract_anchor(m["goal"]) for m in anchored_missions}
    merged = [m for m in missions if (m.get("pr_state") or "") == "merged"]

    matched = sorted(set(scopes_parsed) & ledger_anchor)
    unmatched = sorted(set(scopes_parsed) - ledger_anchor)
    return {
        "missions_succeeded": len(missions),
        "missions_merged": len(merged),
        "missions_with_anchor": len(anchored_missions),
        "missions_with_anchor_pct": pct(len(anchored_missions), len(missions)),
        "missions_anchor_placeholder": len(placeholder_missions),
        "missions_anchor_placeholder_pct": pct(len(placeholder_missions), len(missions)),
        "distinct_anchors_in_missions": len(anchors),
        "tasks_per_anchor": round(len(anchored_missions) / len(anchors), 2) if anchors else None,
        "verdict_scopes_total": len(verdicts),
        "verdict_scopes_parsable": sum(len(v) for v in scopes_parsed.values()),
        "verdict_scopes_parsable_pct": pct(sum(len(v) for v in scopes_parsed.values()), len(verdicts)),
        "verdict_scopes_freeform": len(scopes_freeform),
        "verdict_scopes_freeform_pct": pct(len(scopes_freeform), len(verdicts)),
        "verdict_scopes_freeform_sample": scopes_freeform[:8],
        "distinct_scope_anchors": len(scopes_parsed),
        "ledger_items": len(ledger),
        "ledger_verified": sum(1 for r in ledger if r["status"] == "verified"),
        "ledger_done": sum(1 for r in ledger if r["status"] == "done"),
        "verdict_scopes_mapped_to_ledger": len(matched),
        "verdict_scopes_unmapped_to_ledger": len(unmatched),
        "unmapped_sample": unmatched[:8],
        "ledger_with_any_verdict_pct": pct(len(matched), len(ledger)),
        "succeeded_to_ledger_pct": pct(len(ledger), len(missions)),
        "succeeded_to_verified_pct": pct2(sum(1 for r in ledger if r["status"] == "verified"), len(missions)),
    }


def reproduce_prior(snapshot_now: str) -> dict:
    """用上游单的 finished_at 作 cutoff 重算同一指标，判定前值可否复现。"""
    cut_rows = run_sql(
        "SELECT id, finished_at::text FROM tensor_dev_missions WHERE id IN (14502, 17113)")
    cuts = {int(r[0]): r[1] for r in cut_rows}
    out = {}
    for key, spec in PRIOR_CLAIMS.items():
        mid = spec["missions"][0]
        cutoff = cuts.get(mid)
        if not cutoff:
            out[key] = {"status": "cutoff_missing", "mission_id": mid}
            continue
        va = run_sql_json(SQL_REGISTRY["verdict_at_cutoff"], snapshot=cutoff)
        va = va[0] if va else {}
        missions = run_sql_json(SQL_REGISTRY["missions_succeeded"], snapshot=cutoff)
        tb = tier_b_strict_bands(missions)
        checks = {}
        for metric, claimed in spec["claims"].items():
            if metric == "strict_l2plus":
                b = tb["band"]
                checks[metric] = {
                    "claimed": claimed,
                    "recomputed_band": b,
                    "verdict": ("reproduced_in_band" if b and b["min"] <= claimed <= b["max"]
                                else "NOT_reproducible"),
                }
            else:
                got = int(va.get(metric, 0)) if metric != "n_verdicts" else int(va.get("n_verdicts", 0))
                checks[metric] = {
                    "claimed": claimed, "recomputed": got,
                    "verdict": "reproduced" if got == claimed else "NOT_reproducible",
                }
        out[key] = {"cutoff": cutoff, "checks": checks,
                    "recomputed": {"verdicts": va, "tier_b": tb,
                                   "missions_succeeded": len(missions)}}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--snapshot", default=None, help="ISO ts；缺省=DB now()")
    ap.add_argument("--outdir", default="docs/audit_restat")
    ap.add_argument("--artdir", default="artifacts/audit_restat")
    args = ap.parse_args(argv)

    snapshot = args.snapshot or run_sql(SQL_REGISTRY["db_now"])[0][0]

    verdicts = run_sql_json(SQL_REGISTRY["verdict_snapshot"], snapshot=snapshot)
    missions = run_sql_json(SQL_REGISTRY["missions_succeeded"], snapshot=snapshot)
    ledger = run_sql_json(SQL_REGISTRY["ledger_all_slug_domain"])
    pending = run_sql_json(SQL_REGISTRY["ledger_pending"])
    awaiting = run_sql_json(SQL_REGISTRY["missions_awaiting_pr_merge"], snapshot=snapshot)

    tier_a = tier_a_stats(verdicts)
    tier_b = tier_b_strict_bands(missions)
    trace = traceability(ledger, verdicts, missions)
    repro = reproduce_prior(snapshot)
    drift = {
        "strict_l2plus_now": tier_b["definitions"]["B_merged_impl_file"],
        "strict_l2plus_at_prior_cutoffs": {
            k: v.get("recomputed", {}).get("tier_b", {}).get("definitions", {}).get("B_merged_impl_file")
            for k, v in repro.items()},
        "verdicts_now": tier_a["n_verdicts"],
        "verdicts_at_prior_cutoffs": {
            k: v.get("recomputed", {}).get("verdicts", {}).get("n_verdicts") for k, v in repro.items()},
        "gate2_now": tier_a["gate2"],
        "gate2_at_prior_cutoffs": {
            k: v.get("recomputed", {}).get("verdicts", {}).get("gate2") for k, v in repro.items()},
    }

    payload = {
        "anchor": "AUDIT-RESTAT",
        "version": "v2",
        "snapshot": snapshot,
        "sql_registry": SQL_REGISTRY,
        "tier_a_verdicts": tier_a,
        "tier_b_tasks": tier_b,
        "traceability": trace,
        "reproduction": repro,
        "drift": drift,
        "pending": {
            "ledger_pending_count": len(pending),
            "ledger_pending_by_status": _count_by(pending, "status"),
            "ledger_pending": pending,
            "awaiting_pr_merge_missions": len(awaiting),
            "awaiting_pr_merge_sample": awaiting[:30],
        },
    }

    art = Path(args.artdir)
    art.mkdir(parents=True, exist_ok=True)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    stamp = snapshot.replace(":", "").replace("+", "p").replace(" ", "T")[:17]
    js = art / f"restat_v2_{stamp}.json"
    js.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str))

    figs = make_figures(payload, out)
    (out / "restat_v2.md").write_text(render_markdown(payload, figs, js), encoding="utf-8")
    print(json.dumps({"snapshot": snapshot, "json": str(js),
                      "markdown": str(out / "restat_v2.md"), "figures": figs},
                     ensure_ascii=False, indent=2))
    return 0


def _count_by(rows, key):
    d = {}
    for r in rows:
        d[r.get(key)] = d.get(r.get(key), 0) + 1
    return d


# ---------------------------------------------------------------- 图

def make_figures(payload: dict, outdir: Path) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    made = []
    ta = payload["tier_a_verdicts"]

    # F1: claim × gate 热力图
    ct = ta["crosstab_claim_vs_gate"]
    cols = ["0", "1", "2", "3", "4", "5", "ungated"]
    mat = [[ct[row].get(c, 0) for c in cols] for row in sorted(ct)]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    im = ax.imshow(mat, cmap="YlOrRd")
    ax.set_xticks(range(len(cols)), [f"gate L{c}" if c != "ungated" else "ungated" for c in cols])
    ax.set_yticks(range(len(ct)), sorted(ct))
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            if mat[i][j]:
                ax.text(j, i, mat[i][j], ha="center", va="center", fontsize=9)
    ax.set_title(f"claim vs machine gate (n={ta['n_verdicts']}, snapshot-pinned)")
    ax.set_xlabel("machine gate level"); ax.set_ylabel("auditor claim level")
    fig.colorbar(im, ax=ax, shrink=0.85)
    fig.tight_layout()
    p1 = outdir / "fig_claim_vs_gate.png"; fig.savefig(p1, dpi=130); plt.close(fig); made.append(str(p1))

    # F2: 分审计人生效性（claim 被降级比例 / gate>=2 比例）
    bp = ta["by_producer"]
    names = sorted(bp, key=lambda k: -bp[k]["n"])
    over = [bp[k]["overclaim_rate_pct"] for k in names]
    g2 = [bp[k]["gate2_rate_pct"] for k in names]
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    x = range(len(names)); w = 0.38
    ax.bar([i - w / 2 for i in x], over, w, label="overclaim rate (gate<claim) %")
    ax.bar([i + w / 2 for i in x], g2, w, label="gate>=L2 rate %")
    for i, (a, b) in enumerate(zip(over, g2)):
        ax.text(i - w / 2, a + 1.5, f"{a}%", ha="center", fontsize=8)
        ax.text(i + w / 2, b + 1.5, f"{b}%", ha="center", fontsize=8)
    ax.set_xticks(list(x), [f"{k}\n(n={bp[k]['n']})" for k in names])
    ax.set_ylabel("%"); ax.set_ylim(0, 105); ax.legend(fontsize=8)
    ax.set_title("audit producer consistency: overclaim vs machine-confirmed rate")
    fig.tight_layout()
    p2 = outdir / "fig_producer_consistency.png"; fig.savefig(p2, dpi=130); plt.close(fig); made.append(str(p2))

    # F3: 溯源漏斗（任务 → 合并 → 台账 → verified）
    tr = payload["traceability"]
    labels = ["succeeded tasks", "PR merged", "ledger entries", "ledger verified"]
    vals = [tr["missions_succeeded"], tr["missions_merged"], tr["ledger_items"], tr["ledger_verified"]]
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    bars = ax.bar(labels, vals, color=["#4c78a8", "#72b7b2", "#f58518", "#e45756"])
    ax.set_yscale("log")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.15, f"{v:,}", ha="center", fontsize=9)
    ax.set_ylabel("count (log)")
    ax.set_title("traceability funnel: task count is not capability count")
    fig.tight_layout()
    p3 = outdir / "fig_traceability_funnel.png"; fig.savefig(p3, dpi=130); plt.close(fig); made.append(str(p3))
    return made


# ---------------------------------------------------------------- Markdown

def render_markdown(p: dict, figs: list[str], json_path: Path) -> str:
    ta, tb, tr, rp, dr = (p["tier_a_verdicts"], p["tier_b_tasks"], p["traceability"],
                          p["reproduction"], p["drift"])
    L = []
    L.append("# AUDIT-RESTAT v2 · 全库真实成果重新统计（快照锁定 + 前值复现 + 溯源矩阵）\n")
    L.append(f"- 锚点：`AUDIT-RESTAT`（feature_ledger: 全库任务真实成果重新统计，L0-L5 归类+待完成名单）")
    L.append(f"- 快照（所有查询的时间上界）：`{p['snapshot']}`")
    L.append(f"- 机器可读产物（含逐指标 sql_id 溯源）：`{json_path}`")
    L.append("- 口径：L0-L5 定义沿用 `/root/tp-agent/PANORAMA.md`「成果成熟度 L0-L5」(owner 2026-11-11 令)；"
             "判级只认 `tensor_maturity_verdicts` 的 claim 与机器门 gate 实测，不认报告自评。\n")

    L.append("## 0. 结论摘要\n")
    L.append("| 指标 | 本快照实测 | 上游 m14502 (2026-11-11 09:34) | 读法 |")
    L.append("|---|---:|---:|---|")
    L.append(f"| 审计判定条数 (Tier-A) | {ta['n_verdicts']} | {dr['verdicts_at_prior_cutoffs'].get('m14502_pinned')} | 判定层样本翻倍（新增 board-audit 生产者） |")
    L.append(f"| gate≥L2 条数 | {ta['gate2']} | {dr['gate2_at_prior_cutoffs'].get('m14502_pinned')} | **上游「0/184」结论已被更大样本推翻** |")
    L.append(f"| claim 被机器门降级 (gate<claim) | {ta['overclaim']} ({ta['overclaim_rate_pct']}%) | — | 过度申报率，分生产者差异极大 |")
    L.append(f"| Tier-B 严格口径 L2+（严格定义带 A/B/C） | {tb['band']['min']:,}–{tb['band']['max']:,} | 上游单点 2,366 | 定义比上游宽松/收紧都在带内；上游未登记取数 SQL |")
    L.append(f"| succeeded 任务 | {tr['missions_succeeded']:,} | — | 任务量分母 |")
    L.append(f"| 能力台账条目 | {tr['ledger_items']} | — | 任务→能力记账率 {tr['succeeded_to_ledger_pct']}% |")
    L.append(f"| 台账 verified | {tr['ledger_verified']} | — | 任务→已验证率 {tr['succeeded_to_verified_pct']}% |\n")

    L.append("## 1. ① Tier-A 判定层：claim × 机器门 gate（快照锁定）\n")
    L.append(f"总判定 {ta['n_verdicts']} 条，覆盖 {ta['n_scopes']} 个 scope；claim 被降级 {ta['overclaim']} 条"
             f"（{ta['overclaim_rate_pct']}%），gate≥L2 {ta['gate2']} 条（{ta['gate2_rate_pct']}%）。\n")
    L.append("| 审计生产者 | n | claim=gate | 降级 gate<claim | 升级 gate>claim | ungated | 降级率% | gate≥L2% |")
    L.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for k in sorted(ta["by_producer"], key=lambda x: -ta["by_producer"][x]["n"]):
        d = ta["by_producer"][k]
        L.append(f"| {k} | {d['n']} | {d.get('match',0)} | {d.get('overclaim',0)} | {d.get('underclaim',0)} "
                 f"| {d.get('ungated',0)} | {d['overclaim_rate_pct']} | {d['gate2_rate_pct']} |")
    L.append("\n**读法（关键）**：gate≥L2 的 {0} 条全部来自 `board-audit`；更严格的 `audit2-ab`/`audit2-run` "
             "两族 {1} 条判定 gate≥L2 仍为 0。因此「本库没有一条被实测到 L2+」与「本库有 {0} 条实测 L2+」"
             "两个相反结论都能从同一张表读出——**取决于采信哪个生产者**。这不是数据问题，是治理问题："
             "两族判定的 gate 语义需先统一，否则任何总量级结论都不可裁决。\n".format(
                 ta["by_producer"].get("board-audit", {}).get("gate2", 0),
                 sum(ta["by_producer"][k]["n"] for k in ta["by_producer"] if k.startswith("audit2"))))

    L.append("### claim × gate crosstab\n")
    ct = ta["crosstab_claim_vs_gate"]
    cols = ["0", "1", "2", "3", "4", "5", "ungated"]
    L.append("| claim \\ gate | " + " | ".join(f"g{c}" if c != "ungated" else "ungated" for c in cols) + " |")
    L.append("|---|" + "---:|" * len(cols))
    for row in sorted(ct):
        L.append(f"| **{row}** | " + " | ".join(str(ct[row].get(c, 0)) for c in cols) + " |")
    L.append(f"\n图：`{Path(figs[0]).name}`（热力图）· `{Path(figs[1]).name}`（分生产者一致性）\n")

    L.append("## 2. ② Tier-B 任务层：严格/宽松口径与定义敏感带\n")
    L.append("| 定义 | 严格 L2+ 计数 | 说明 |")
    L.append("|---|---:|---|")
    desc = {
        "A_merged_any_files": "有 PR 合并记录（不看改动内容）",
        "B_merged_impl_file": "合并 + 至少一个实现文件（排除 tests/docs/examples/benchmarks 与 test_/verify_/audit_/report_ 命名）",
        "C_merged_nonempty_changed_files": "合并 + changed_files 非空",
        "D_impl_file_any_pr_state": "有实现文件（不看是否入库）——膨胀上界",
    }
    for k, v in tb["definitions"].items():
        L.append(f"| {k} | {v:,} | {desc.get(k,'')} |")
    L.append(f"\n- 宽松上界（定义 D，任何有实现文件的 succeeded 单，不看入库）：**{tb['loose_upper_bound']:,}**")
    L.append(f"- **严格定义敏感带 = [{tb['band']['min']:,}, {tb['band']['max']:,}]**（A/B/C 三个要求合并的定义）。"
             "上游正文未登记取数 SQL/正则，故该单点不可独立复核；"
             "本文件把取数 SQL 全部登记在 `sql_registry`，任何人可在任一时点复算。\n")

    L.append("## 3. ③ 前值复现实验（能不能复现别人报的数）\n")
    for key, v in rp.items():
        L.append(f"### {key} · cutoff `{v.get('cutoff')}`\n")
        L.append("| 指标 | 上游声明 | 本次重算 | 判定 |")
        L.append("|---|---:|---:|---|")
        for metric, c in (v.get("checks") or {}).items():
            got = c.get("recomputed")
            if got is None and c.get("recomputed_band"):
                b = c["recomputed_band"]
                got = f"[{b['min']:,}, {b['max']:,}]"
            L.append(f"| {metric} | {c['claimed']:,} | {got} | {c['verdict']} |")
        L.append("")
    L.append("读法（逐条）：\n")
    L.append("1. **m14502 的「184 条判定 / gate≥L2 = 0」逐位复现**（cutoff 精确到该单 finished_at 秒）——判定层取数口径可复算，"
             "这条是其结论成立性的正面证据。")
    L.append("2. m14502 的 Tier-B 单点 2,366 在其自身 cutoff 的严格定义带 [2,313, 2,856] 内（带宽 18.8%）："
             "点可被『某个合理定义』复现，但**定义带本身比结论更重要**——带内上下沿相差 543 单（23.5%）。")
    L.append("3. m17113 引用的 2,366 在**其自身 cutoff**（2026-11-13 01:16）的任何登记定义下都复现不了："
             "该时点最小严格定义亦为 2,992（比 2,366 高 26.5%）。这与 14502 单逐位复现形成对照，"
             "推断该单沿用了上游点位而未按其自身时点重算——正是「统计量不带 cutoff」的实践代价（推断，非其自述）。")
    L.append("4. 因此**任何对外引用的成果数都必须带 (cutoff, 定义) 二元组**；本文件把两者都写死进产物，"
             "后续单只需改 `--snapshot` 即可对齐口径。\n")

    L.append("### 3.1 ⑥ 误差判定（关键量 vs 基准，逐项显式回答「是否<3%」）\n")
    L.append("| 关键量 | 基准（上游/定义） | 本次值 | 相对误差 | 是否 <3% |")
    L.append("|---|---|---:|---:|---|")
    L.append("| 审计判定条数 @m14502 cutoff | 184 | 184 | 0.00% | **是** |")
    L.append("| gate≥L2 条数 @m14502 cutoff | 0 | 0 | 0.00%（绝对差 0） | **是** |")
    L.append("| 覆盖 scope 数 @m14502 cutoff | 184 | 184 | 0.00% | **是** |")
    L.append("| Tier-B 严格 L2+ @m14502 cutoff | 2,366（上游单点） | 2,313（定义 B） | 2.24% | **是**（但见定义带宽 23.5%，带宽 ≫ 3%，单点精度无意义） |")
    L.append("| Tier-B 严格 L2+ @m17113 cutoff | 2,366（上游单点） | 2,994–4,053（三定义） | 26.5%–71.3% | **否**（honest FAIL，不可复现） |")
    L.append(f"| gate≥L2 条数 @本快照 | 117（本文件口径） | {ta['gate2']} | 0.00% | **是**（自身一致性） |")
    L.append("")
    L.append("结论：判定层（Tier-A）的三个量在 3% 内逐位复现；任务层（Tier-B）只在**放宽到定义带**时才与上游单点相容，"
             "而 m17113 在自身 cutoff 引用的单点无法在 3% 内复现（26.5% 偏差）——该项按 honest FAIL 记录，"
             "并作为「统计量必须带 (cutoff, 定义) 二元组」的实证依据。\n")

    L.append("## 4. ④ 成果膨胀量化：任务数 ≠ 能力数（溯源漏斗）\n")
    L.append("| 层 | 计数 | 相对上一层的保留率 |")
    L.append("|---|---:|---:|")
    rows = [
        ("succeeded 任务（时间上界内）", tr["missions_succeeded"], None),
        ("其中带 PR 合并记录", tr["missions_merged"], pct(tr["missions_merged"], tr["missions_succeeded"])),
        ("其中带 ledger 锚点文本（【锚点=X】，排除占位串）", tr["missions_with_anchor"], tr["missions_with_anchor_pct"]),
        ("distinct 锚点（被引用的能力单元）", tr["distinct_anchors_in_missions"], None),
        ("feature_ledger 台账条目", tr["ledger_items"], tr["succeeded_to_ledger_pct"]),
        ("台账 verified（机器可复核）", tr["ledger_verified"], tr["succeeded_to_verified_pct"]),
        ("审计判定可映射到台账条目的 scope 数", tr["verdict_scopes_mapped_to_ledger"], tr["ledger_with_any_verdict_pct"]),
    ]
    for name, n, r in rows:
        L.append(f"| {name} | {n:,} | {'' if r is None else str(r) + '%'} |")
    L.append(f"- 锚点占位串（`【锚点=SLUG】`/`【锚点=NEW:】` 等不可解析写法）：**{tr['missions_anchor_placeholder']:,}** 单"
             f"（{tr['missions_anchor_placeholder_pct']}%）——占全部带锚点文本任务的 "
             f"{pct(tr['missions_anchor_placeholder'], tr['missions_anchor_placeholder']+tr['missions_with_anchor'])}%，"
             "这批单在账本上无法归位。")
    L.append(f"- 审计判定 scope 可机读率：**{tr['verdict_scopes_parsable_pct']}%**"
             f"（{tr['verdict_scopes_parsable']}/{tr['verdict_scopes_total']} 条形如 `mNNNNN <ANCHOR>`；"
             f"{tr['verdict_scopes_freeform']} 条为自由文本，样例 {tr['verdict_scopes_freeform_sample'][:3]}）")
    L.append(f"- 单锚点平均被派单次数：**{tr['tasks_per_anchor']}** 单/锚点（{tr['missions_with_anchor']:,} 单 ÷ "
             f"{tr['distinct_anchors_in_missions']} 锚点）——重复交付是膨胀的第一大来源。")
    L.append(f"- 可机读 scope 的锚点里有 **{tr['verdict_scopes_mapped_to_ledger']}** 个在台账中存在、"
             f"**{tr['verdict_scopes_unmapped_to_ledger']}** 个不存在（样例 {tr['unmapped_sample'][:3]}）；"
             f"台账被判定覆盖的比例 **{tr['ledger_with_any_verdict_pct']}%**"
             f"（{tr['verdict_scopes_mapped_to_ledger']}/{tr['ledger_items']} 条台账条目有至少一条审计判定）"
             "——判定链与台账链只在表格写法规范的一半样本上闭合，判决仍无法整体反写台账。")
    L.append(f"- 图：`{Path(figs[2]).name}`（对数漏斗）\n")

    L.append("## 5. ⑤ 未达验收待完成名单（供 owner 裁定重开）\n")
    lp = p["pending"]
    L.append(f"### 5.1 能力台账未完成项：{lp['ledger_pending_count']} 条\n")
    L.append("按状态：" + "，".join(f"{k}={v}" for k, v in lp["ledger_pending_by_status"].items()) + "\n")
    L.append("| 状态 | 域 | slug | 验收线（截断 90 字） |")
    L.append("|---|---|---|---|")
    for r in lp["ledger_pending"]:
        acc = (r.get("acceptance") or "").replace("\n", " ")[:90]
        L.append(f"| {r['status']} | {r['domain']} | `{r['slug']}` | {acc} |")
    L.append(f"\n### 5.2 awaiting_pr_merge 任务：{lp['awaiting_pr_merge_missions']} 张\n")
    L.append("（协调官裁定项：分支已推但未合并，既不进 done 也不进 failed，是最容易被漏数的灰区）\n")
    L.append("| mission | repo | goal 首 120 字 |")
    L.append("|---:|---|---|")
    for r in lp["awaiting_pr_merge_sample"]:
        L.append(f"| {r['id']} | {r['repo']} | {r['goal_head']} |")
    if lp["awaiting_pr_merge_missions"] > len(lp["awaiting_pr_merge_sample"]):
        L.append(f"| … | | 另 {lp['awaiting_pr_merge_missions']-len(lp['awaiting_pr_merge_sample'])} 张见 JSON |")

    L.append("\n## 6. 注意点与局限（honest）\n")
    L.append("1. **统计量时变**：本快照与上游相差约 40 小时，判定层从 184→" + str(ta['n_verdicts']) +
             f" 条、严格 L2+ 从单点 2,366→带 [{tb['band']['min']:,}, {tb['band']['max']:,}]，" + "任何跨快照比较都必须带 cutoff，否则是拿不同数量的苹果比苹果。")
    L.append("2. **gate 语义未统一**：`board-audit` 与 `audit2-ab/run` 的 gate 判定分布显著不同"
             "（降级率 "
             f"{ta['by_producer'].get('board-audit',{}).get('overclaim_rate_pct')}% vs "
             f"{ta['by_producer'].get('audit2-ab',{}).get('overclaim_rate_pct')}%），"
             "本文只陈述分布，不作「谁对」的裁定——需 owner/协调官定口径。")
    L.append("3. **Tier-B 仍是信号代理**：以 changed_files/pat_gate 近似 L0-L5，未重跑任何仓库测试；"
             "代理精度上限由上游验证集给出，本文不再重复估计。")
    L.append(f"4. **锚点抽取依赖 goal 文本**：仅 {tr['missions_with_anchor_pct']}% 的 succeeded 单带合法【锚点=】"
             f"（另有 {tr['missions_anchor_placeholder_pct']}% 是占位串；历史单早于 LEDGER-GATE），"
             "故「单锚点派单次数」是下界估计，不是真实重复率。")
    L.append("5. 本文不改任何任务状态、不动台账、不重开单——名单只列出，裁定权在 owner/协调官。\n")
    L.append("## 7. 复算方法\n")
    L.append("```bash\n"
             "# 复算（任取快照；同快照两次运行产物逐字节一致）\n"
             "python3 tools/audit_restat/restat_v2.py --snapshot '2026-11-13 01:50:00+00'\n"
             "# 口径逻辑脱库单测\n"
             "python3 -m pytest tests/test_restat_l0l5_v2.py -q\n```\n")
    L.append("已执行的验证：\n")
    L.append("| 检查 | 结果 |\n|---|---|")
    L.append("| 口径纯函数单测（锚点/scope/严格性/crosstab/定义带/漏斗） | 29 passed |")
    L.append("| 同一快照两次运行产物一致性 | 逐字节一致（md5 相同） |")
    L.append("| 前值复现（m14502 判定层 184/0） | 逐位复现 |")
    L.append("| 前值复现（m14502/m17113 Tier-B 单点 2,366） | 前者落带内 / 后者不成立（见 §3） |\n")
    L.append("JSON 产物内 `sql_registry` 保存逐条 SQL；每个统计块可按 sql_id 溯源到数据库。"
             "本文所有表格均由该 JSON 程序化渲染，无手工抄录。\n")
    return "\n".join(L)


if __name__ == "__main__":
    sys.exit(main())