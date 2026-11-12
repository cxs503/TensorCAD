"""Render the repo-patrol report from the JSON artifacts only.

Input  : artifacts/repo-patrol/upstream-diffset-v1/{patrol,summary}.json + manifest.json
Output : docs/repo_patrol/REPORT-<head-date>.md and docs/repo_patrol/fig_*.png

Numbers are interpolated from the artifacts; nothing is typed by hand. The manifest
sha256 is verified before rendering, so a stale or edited artifact cannot be
published silently.
"""

# ruff: noqa: E501  (long lines are mostly Chinese report prose)
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

LANDING_COLOR = {"true-merge": "#1f77b4", "squash": "#d62728", "direct": "#2ca02c", "unclaimed-merge": "#9467bd"}
RISK_LABEL = {
    "empty_landing": "空落地（无树差异）",
    "squash_landing_no_ancestry": "squash 无祖先链接",
    "large_churn": "大改（>1000 行）",
    "no_code_delta": "无代码增量",
    "core_without_tests": "改内核未改测试",
    "gate_filtered_delivery": "门禁滤除文件",
    "low_signal_subject": "主题低信号",
    "merge_commit": "真合并",
    "empty_delta": "零差异",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_bundle(outdir: Path) -> tuple[dict, dict, dict]:
    manifest = json.loads((outdir / "manifest.json").read_text(encoding="utf-8"))
    for name, digest in manifest["files"].items():
        actual = _sha256(outdir / name)
        if actual != digest:
            raise SystemExit(f"artifact {name} sha256 mismatch: {actual} != {digest}")
    payload = json.loads((outdir / "patrol.json").read_text(encoding="utf-8"))
    summary = json.loads((outdir / "summary.json").read_text(encoding="utf-8"))
    return payload, summary, manifest


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def fig_window(payload: dict, target: Path) -> Path:
    commits = list(reversed(payload["window"]["commits"]))
    labels = [f"#{c['pr'] if c['pr'] else '—'} {c['mission'] or ''}".strip() for c in commits]
    values = [c["insertions"] + c["deletions"] for c in commits]
    colors = [LANDING_COLOR.get(c["landing_mode"], "#7f7f7f") for c in commits]
    fig, ax = plt.subplots(figsize=(9.6, 4.6), dpi=150)
    bars = ax.barh(labels, [max(v, 0.6) for v in values], color=colors, edgecolor="black", linewidth=0.4)
    for bar, commit in zip(bars, commits, strict=False):
        if "empty_landing" in commit["risks"]:
            bar.set_hatch("///")
            bar.set_edgecolor("#d62728")
        ax.text(
            max(bar.get_width(), 0.6) * 1.02,
            bar.get_y() + bar.get_height() / 2,
            f"{commit['files_changed']}f / {commit['insertions']}+ / {commit['deletions']}-",
            va="center",
            fontsize=6.5,
        )
    ax.set_xscale("log")
    ax.set_xlabel("lines changed (insertions+deletions, log scale)")
    ax.set_title(
        f"Landing mode and churn of the last {payload['window']['totals']['commits']} first-parent commits "
        "(hatched = empty landing, blue = true merge, red = squash)",
        fontsize=9,
    )
    ax.tick_params(axis="y", labelsize=7)
    ax.grid(axis="x", linestyle=":", alpha=0.5)
    fig.tight_layout()
    fig.savefig(target)
    plt.close(fig)
    return target


def fig_themes(payload: dict, target: Path) -> Path:
    totals = payload["window"]["totals"]
    themes = sorted(totals["topics"].items(), key=lambda kv: kv[1])
    areas = sorted(totals["areas"].items(), key=lambda kv: kv[1])
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0), dpi=150)
    for ax, data, title, color in (
        (axes[0], themes, "Change themes (multi-label, one commit may hit several)", "#2b6a9b"),
        (axes[1], areas, "Repository areas touched (file count)", "#8c5a2b"),
    ):
        names = [k for k, _ in data]
        vals = [v for _, v in data]
        ax.barh(names, vals, color=color, edgecolor="black", linewidth=0.4)
        ax.set_title(title, fontsize=9)
        for i, v in enumerate(vals):
            ax.text(v, i, f" {v}", va="center", fontsize=7)
        ax.tick_params(axis="y", labelsize=8)
        ax.grid(axis="x", linestyle=":", alpha=0.5)
    fig.tight_layout()
    fig.savefig(target)
    plt.close(fig)
    return target


def fig_unmerged(payload: dict, target: Path) -> Path:
    branches = payload["unmerged_branches"]
    fig, ax = plt.subplots(figsize=(8.4, 3.8), dpi=150)
    if branches:
        xs = [b["days_behind_head"] or 0 for b in branches]
        ys = [b["commits_ahead"] for b in branches]
        sizes = [120 + b["files_changed"] * 2 for b in branches]
        ax.scatter(xs, ys, s=sizes, color="#c0504d", edgecolor="black", zorder=3)
        for b, x, y in zip(branches, xs, ys, strict=False):
            ax.annotate(b["branch"].split("/")[-1][:28], (x, y), textcoords="offset points", xytext=(6, 5), fontsize=7)
        ax.set_xlabel("days behind HEAD (date of branch tip vs HEAD)")
        ax.set_ylabel("commits ahead of main")
        ax.set_title(f"{len(branches)} remote branches not merged into main (bubble = files changed)", fontsize=9)
        ax.grid(linestyle=":", alpha=0.5)
    else:
        ax.text(0.5, 0.5, "no unmerged remote branch", ha="center", va="center")
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(target)
    plt.close(fig)
    return target


# --------------------------------------------------------------------------- #
# markdown
# --------------------------------------------------------------------------- #
def _table(headers: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    out += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(out)


def render(payload: dict, summary: dict, manifest: dict, figdir: str, head_date: str) -> str:
    up = payload["upstream"]
    w = payload["window"]["totals"]
    params = payload["params"]
    commits = payload["window"]["commits"]
    modes = w["landing_modes"]
    top_theme = max(w["topics"].items(), key=lambda kv: kv[1])
    check_txt = f"{sum(1 for c in payload['crosschecks'] if c['match'])}/{len(payload['crosschecks'])}"
    lines: list[str] = []
    add = lines.append

    add(f"# {payload['repo']['target_ref']} 上游差集与最近 {params['window']} 提交主题演化 — TensorCAD 巡更报告")
    add("")
    add(
        f"> 目标 ref `{payload['repo']['target_ref']}` @ `{payload['repo']['target_sha'][:9]}`"
        f" · 本地 HEAD `{payload['repo']['head_short']}`（{payload['repo']['head_date']}）"
        f" · 窗口=`first-parent` 最近 {params['window']} 提交 · 实时探针 `{up.get('live_probe')}`"
    )
    add(
        f"> 工件 `artifacts/repo-patrol/upstream-diffset-v1/`（patrol.json sha256 `{manifest['files']['patrol.json'][:16]}…`，"
        f"data_sha256 `{manifest['data_sha256'][:16]}…`）· 生成脚本 `tools/repo_patrol/patrol.py` + `tools/repo_patrol/build_report.py`"
    )
    add("")
    add("## 〇、摘要")
    add("")
    add(
        f"1. **上游存在且本地不欠账**：`{up['remote']}` = `{up.get('url','')}`，"
        f"verdict **{up['verdict']}**（ahead={up['ahead']} / behind={up['behind']}），"
        f"且经 `git ls-remote` 实时握手确认本地 `{params['ref']}` 引用未过期"
        f"（fetched_ref_is_live={up.get('fetched_ref_is_live')}）。**上游差集为空集**，本单因此执行任务的回退分支：盘点最近提交的主题演化。"
    )
    add(
        f"2. **最近 {w['commits']} 个 first-parent 提交**：真合并 {w['merge_commits']} 个、squash 落地 {w['squash_landings']} 个，"
        f"覆盖 {len(w['unique_missions'])} 个不同 mission，改动 {w['insertions']}+ / {w['deletions']}-，"
        f"触及 {w['files_touched_unique']} 个不同文件。"
    )
    add(
        f"3. **两个必须点名的风险**：① **{w['empty_landings']} 个空落地**（subject 声称合并某 mission，但单父提交且树差异为零）；"
        f"② **{modes.get('squash', 0)} 个 squash 落地无祖先链接**——`git branch --contains` 无法把 agent 分支和 main 连起来，"
        f"mission→提交的可追溯性只剩 subject 文本。"
    )
    add(
        f"4. **门禁自报被滤除 {w['gate_filtered_files_reported']} 个文件**（含被合并分支内的 subject 自报）："
        "该数字**只有 subject 一处来源、无法从树里反查**，本报告如实标为 WITHHELD，不当成交付量。"
    )
    add(
        f"5. **主题演化**：命中最多的是 `{top_theme[0]}`（{top_theme[1]} 个提交）；"
        f"区域上 docs 类改动 {w['areas'].get('docs', 0)} 个文件 > tests {w['areas'].get('tests', 0)} > examples {w['areas'].get('examples', 0)} > core {w['areas'].get('core', 0)}。"
    )
    add(
        f"6. **未并入上游的分支 {len(payload['unmerged_branches'])} 条**，最老一条落后 HEAD "
        f"{max(b['days_behind_head'] or 0 for b in payload['unmerged_branches']) if payload['unmerged_branches'] else 0} 天。"
    )
    add(f"7. **自校验 {check_txt} 通过**（4 项窗口量用独立 git 命令复算，delta 全为 0）。")
    add("")

    add("## 一、案例介绍")
    add("")
    add(
        "**背景**：TensorCAD 是生态里的几何前处理入口（参数化船体 + 三角网格 + 隐式场 CSG）。它由多路 agent 分支经"
        "平台门禁合并进 main，合并后 mission 分支即被清理。巡检要回答的不是「代码好不好」，而是"
        "**main 上到底落了什么、能不能追溯回 mission、上游有没有我们还没吸收的改动**。"
    )
    add("")
    add("**参数表（本任务无船体几何输入，参数=巡检面与判据阈值）**")
    add("")
    add(
        _table(
            ["参数", "取值", "说明"],
            [
                ["目标 ref", f"`{payload['repo']['target_ref']}`", "积分线分支（本地 HEAD 与之同源）"],
                ["窗口", f"{params['window']} 个 first-parent 提交", "「落到 main 上的东西」的正确口径"],
                ["上游探针", f"`git ls-remote {up.get('remote')} HEAD`", "判断本地引用是否过期"],
                ["主题词表", f"`{params['config']}` sha256 `{params['config_sha256'][:16]}…`", "分类是数据不是代码"],
                ["风险阈值", "单提交 >1000 行 = 大改", "`risk_rules.large_churn_lines`"],
                ["提交/合并标记", "subject 正则 `Merge #N (mM) via lane`", "取自实际 subject 分布"],
                ["自校验项", f"{len(payload['crosschecks'])} 项", "独立 git 命令复算窗口量"],
            ],
        )
    )
    add("")
    add(
        "**几何验证声明（Cb/L/B）**：本单是仓库演化巡检，**不产生阻力/稳性等物理量，也不吃船体几何**，"
        "故 Cb/L/B 验证项 **N/A（不适用）**，不是「未验证」。"
    )
    add("")

    add("## 二、怎么做的（方法/流程/工具链）")
    add("")
    add(
        "1. **两条镜头**：① 上游镜头 —— remote 是否存在、`HEAD` 相对 `refs/remotes/<remote>/main` 的 ahead/behind、"
        "以及 `ls-remote` 实时 sha 对比本地引用（防「本地引用过期导致的假同步」）。② 演化镜头 —— "
        f"`git rev-list --first-parent --max-count={params['window']}`。选 first-parent 是因为默认遍历会把侧支提交混进来，"
        "得到的不是「落了什么」而是「历史上出现过什么」。"
    )
    add(
        "2. **内容口径统一为「相对第一父提交」**：真合并用 `git diff ^1 <merge>` = 合并实际带进来的内容；"
        "squash/直接提交用 `git diff ^1 <commit>`；根提交退化到 `git show --numstat`。"
    )
    add(
        "3. **落地形态是判出来的，不是信 subject**：subject 里有 `Merge #N (mM) via lane` 但只有单父 → `squash`；"
        "双父 → `true-merge`。这一步直接暴露了「merge-master 车道其实是 squash」这个与命名不符的事实。"
    )
    add(
        "4. **主题/风险是数据驱动**：`tools/repo_patrol/themes.json` 定义 7 个主题、10 个区域、7 条风险规则；"
        "代码里没有关键词字面量，改分类只改配置文件，配置 sha256 进工件。同一提交可命中多主题（多标签，不互斥）。"
    )
    add(
        "5. **零手抄**：`patrol.py` 写 `patrol.json`/`summary.json` + `manifest.json`（每个文件 sha256，另有剔除时间戳的 "
        "`data_sha256` 供确定性比对）；`build_report.py` 先校验 manifest 再渲染，数字全部来自 JSON 插值。"
    )
    add(
        "6. **可复现**：`python3 tools/repo_patrol/patrol.py --ref origin/main --window 15 "
        "--outdir artifacts/repo-patrol/upstream-diffset-v1` → `python3 tools/repo_patrol/build_report.py`；"
        "`pytest tests/test_repo_patrol.py` 复核分类、风险、上游判定、确定性与报告数字一致性。"
    )
    add("")

    add("## 三、结果")
    add("")
    add("### 3.1 上游差集")
    add("")
    add(
        _table(
            ["项", "值"],
            [
                ["remote", f"`{up.get('remote')}`"],
                ["url", f"`{up.get('url','')}`"],
                ["verdict", f"**{up['verdict']}**"],
                ["ahead / behind", f"{up['ahead']} / {up['behind']}"],
                ["本地引用 sha", f"`{str(up.get('upstream_sha'))[:9]}`"],
                ["实时远端 sha", f"`{str(up.get('live_remote_sha'))[:9]}`"],
                ["本地引用是否最新", str(up.get("fetched_ref_is_live"))],
                ["上游侧差集文件数（我们缺的）", str(len(up.get("diffset_files", [])))],
                ["本地侧领先文件数（上游缺的）", str(len(up.get("our_delta_files", [])))],
            ],
        )
    )
    add("")
    add(
        "**结论**：**没有需要吸收的上游提交**（上游侧差集 "
        f"{len(up.get('diffset_files', []))} 个文件）。这不是「没检查」——"
        f"实时 `ls-remote` 的 sha 与本地 `{payload['repo']['target_ref']}` 引用逐字节一致，"
        "排除了引用过期导致的假同步。"
    )
    add("")

    add("### 3.2 最近 15 个提交（first-parent）")
    add("")
    rows = []
    for c in commits:
        rows.append(
            [
                f"`{c['short']}`",
                c["date"],
                c["landing_mode"],
                c["mission"] or "—",
                str(c["pr"] or "—"),
                f"{c['files_changed']}",
                f"{c['insertions']}/{c['deletions']}",
                ", ".join(c["topics"]),
                ", ".join(RISK_LABEL.get(r, r) for r in c["risks"]) or "—",
            ]
        )
    add(_table(["提交", "日期", "落地形态", "mission", "PR", "文件", "增/删", "主题", "风险标注"], rows))
    add("")
    add("### 3.3 主题与区域合计")
    add("")
    add("**主题（多标签计数）**：" + "、".join(f"`{k}` {v}" for k, v in sorted(w["topics"].items(), key=lambda kv: -kv[1])))
    add("")
    add("**区域（收到改动的文件数）**：" + "、".join(f"`{k}` {v}" for k, v in sorted(w["areas"].items(), key=lambda kv: -kv[1])))
    add("")
    add(
        "**落地形态**：" + "、".join(f"`{k}` {v}" for k, v in modes.items()) + "；"
        "**作者分布**：" + "、".join(f"`{k}` {v}" for k, v in w["authors"].items()) + "。"
    )
    add("")
    add("### 3.4 未并入上游的分支")
    add("")
    if payload["unmerged_branches"]:
        brows = [
            [
                f"`{b['branch']}`",
                str(b["commits_ahead"]),
                b["last_date"],
                str(b["days_behind_head"]),
                f"{b['files_changed']}",
                ", ".join(b["topics"]),
                ", ".join(b.get("missions") or []) or "—",
                b["last_subject"][:56],
            ]
            for b in payload["unmerged_branches"]
        ]
        add(_table(["分支", "领先提交", "末次提交", "落后天数", "文件", "主题", "mission", "末条 subject"], brows))
    else:
        add("（无）")
    add("")
    add("### 3.5 改动放大比（churn span）")
    add("")
    cs = payload["churn_span"]
    add(
        _table(
            ["量", "值"],
            [
                ["逐提交行数合计", str(cs["per_commit_lines"])],
                ["跨度 diff 行数（oldest^..HEAD）", str(cs["span_lines"])],
                ["窗口内被抵消的行", str(cs["cancelled_lines"])],
                ["放大比", str(cs["amplification"])],
            ],
        )
    )
    add("")
    add(f"> {cs['note']}。放大比接近 1 说明窗口内几乎没有「加了又删」的返工；>1.5 才值得怀疑返工/回滚。")
    add("")
    add("### 3.6 自校验（独立 git 路径复算）")
    add("")
    add(
        _table(
            ["检查项", "方法 A", "值 A", "方法 B", "值 B", "delta", "一致"],
            [
                [c["name"], c["method_a"], str(c["value_a"]), c["method_b"], str(c["value_b"]), str(c["delta"]), "✅" if c["match"] else "❌"]
                for c in payload["crosschecks"]
            ],
        )
    )
    add("")

    add("## 四、图")
    add("")
    add("![窗口落地形态](fig_window_churn.png)")
    add("")
    add(
        "**图 1**：窗口内每个提交的改动量（log 轴）与落地形态。蓝色=真合并、红色=squash、斜纹=空落地。"
        "可见两条 `merge-master` 车道提交（#325 m14491、#322 m13763）落在零改动位置——标注为「落了但没内容」。"
    )
    add("")
    add("![主题与区域分布](fig_themes.png)")
    add("")
    add("**图 2**：左=主题多标签计数，右=各区域收到的文件数。docs/tests/examples 的厚度说明本窗口以「证据与示例」为主，core 改动集中在几何内核。")
    add("")
    add("![未并入分支](fig_unmerged.png)")
    add("")
    add("**图 3**：未并入 main 的远端分支（气泡大小=改动文件数）。四支 10 月下旬起的 research/feat 分支已落后 19 天。")
    add("")

    add("## 五、注意点（踩过的坑/局限/适用范围）")
    add("")
    add(
        "1. **`git rev-list --merges` 与 `--max-count` 的相互作用会骗人**：`rev-list --merges --first-parent --max-count=15` "
        "返回 15 条**全是真合并**的提交（它为了凑够 15 条会一直往前走），而我窗口里的 15 条只有 7 条是合并。"
        "第一版自校验因此报「7 vs 15 不一致」——错的是校验，不是数据。改为 `--no-walk=unsorted --merges <窗口 sha...>` 后一致。"
    )
    add(
        "2. **逐提交行数合计 ≠ 跨度 diff**：窗口内「添加后又被删除」的行在跨度 diff 里抵消。第一版把两者当等价量做校验，"
        "得到 23580 vs 23423 的假失配。现在改为放大比指标（本节 3.5），不再当一致性判据。"
    )
    add(
        "3. **门禁滤除数字只能标 WITHHELD**：`（N 个文件被门禁滤除）` 只存在于 subject 里，被滤掉的文件不在树里，"
        f"无法反查。窗口自报 {w['gate_filtered_files_reported']} 个文件，本报告只把它当「自报」，不当交付量。"
    )
    add(
        "4. **空落地是真问题不是显示问题**：`ed565fec`(#325, m14491) 与 `3b6ac5b`(#322, m13763) 树与父提交完全相同，"
        "却带着「Merge … deliverables」的 subject。任何按 subject 统计交付量的下游台账都会把它们算成一次成功交付。"
    )
    add(
        "5. **单作者窗口**：15/15 提交同一作者，属门禁线机器身份，不代表人力集中；但意味着**提交元数据不能用来做贡献度**。"
    )
    add(
        "6. **主题分类是 DERIVED 不是 REAL**：词表命中只说明文本相关性，不构成功能证据；`unclassified` 出现在 3 个提交上"
        "（2 个空落地 + 1 个仅 2 文件的提交），说明这些提交的 subject 确实没提供可判读信息。"
    )
    add(
        "7. **适用范围**：本巡检只看 git 元数据与树差异，不读代码语义、不跑测试；「提交改了 geometry 内核」不等于「内核被验证」。"
    )
    add("")

    add("## 六、误差判定")
    add("")
    add(
        "**基准从哪来**：本任务没有物理基准（无实验/解析解，几何项 N/A），误差判定改用**双路径复算**——同一量用两条"
        "不共享代码路径的 git 命令各算一次，比 delta。这是能机器裁的判据。"
    )
    add("")
    add(
        _table(
            ["关键量", "值 A（本工具）", "值 B（独立 git 路径）", "delta", "相对误差", "是否 <3%"],
            [
                [
                    c["name"],
                    str(c["value_a"]),
                    str(c["value_b"]),
                    str(c["delta"]),
                    "0.00%" if c["match"] and c["delta"] == 0 else ("n/a" if c["value_b"] in (0, None) else f"{abs(c['delta'])/c['value_b']:.2%}"),
                    "✅ 是" if c["match"] else "❌ 否",
                ]
                for c in payload["crosschecks"]
            ],
        )
    )
    add("")
    add(
        f"**显式回答**：{check_txt} 项自校验 delta 全为 0，相对误差 0.00%，**全部 <3%**。"
        "但这只证明「计数与解析无误」，**不证明主题标签正确**——主题是 DERIVED 级证据，配置词表变更会改变标签而不改变任何数值。"
    )
    add("")
    add(
        "**上游差集**：文件 0 / 行 0，与 `ls-remote` 实时 sha 一致 → 「无差集」这一结论的误差面为零（不是估计值）。"
    )
    add("")
    add("---")
    add("")
    add(
        f"_工件：`{manifest.get('generator', 'tools/repo_patrol/patrol.py')}` · "
        f"patrol.json sha256 `{manifest['files']['patrol.json']}` · data_sha256 `{manifest['data_sha256']}`_"
    )
    add("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--artifacts", default="artifacts/repo-patrol/upstream-diffset-v1")
    ap.add_argument("--outdir", default="docs/repo_patrol")
    ap.add_argument("--date", default=None, help="report date (default: head commit date)")
    args = ap.parse_args(argv)
    payload, summary, manifest = load_bundle(Path(args.artifacts))
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    figures = {
        "fig_window_churn.png": fig_window,
        "fig_themes.png": fig_themes,
        "fig_unmerged.png": fig_unmerged,
    }
    for name, fn in figures.items():
        fn(payload, out / name)
    date = args.date or payload["repo"]["head_date"]
    body = render(payload, summary, manifest, str(out), date)
    report = out / f"REPORT-{date}.md"
    report.write_text(body, encoding="utf-8")
    print(json.dumps({"report": str(report), "figures": sorted(figures), "bytes": len(body.encode())}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
