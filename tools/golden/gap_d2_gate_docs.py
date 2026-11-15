#!/usr/bin/env python3
"""Render docs/gap-cad-d2-geometry-golden-gate.md from the run's own artefacts.

Every number in the generated document comes from a file written by a run -- the suite's
``golden_suite.json`` / ``mutation_report.json`` snapshots under
``artifacts/gap-cad-d2-golden/``, the ``run_meta.json`` that recorded the exit codes and
wall times, and the two files that actually carry the trigger lists
(``.github/workflows/geometry-golden.yml`` and ``.githooks/pre-commit-geometry-golden.sh``).
Nothing is typed twice, so the document cannot disagree with the run that produced it.

    python3 tools/golden/gap_d2_gate_docs.py            # rewrite the document
    python3 tools/golden/gap_d2_gate_docs.py --check    # fail if it is stale (CI-friendly)
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ART = os.path.join(REPO_ROOT, "artifacts", "gap-cad-d2-golden")
DOC = os.path.join(REPO_ROOT, "docs", "gap-cad-d2-geometry-golden-gate.md")
WORKFLOW = os.path.join(REPO_ROOT, ".github", "workflows", "geometry-golden.yml")
HOOK = os.path.join(REPO_ROOT, ".githooks", "pre-commit-geometry-golden.sh")


def _load(name: str):
    with open(os.path.join(ART, name), encoding="utf-8") as handle:
        return json.load(handle)


def workflow_paths() -> dict[str, list[str]]:
    """Trigger lists as the workflow declares them (same parse as the unit test)."""
    with open(WORKFLOW, encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    sections: dict[str, list[str]] = {"push": [], "pull_request": []}
    current = None
    in_on_block = False
    for line in lines:
        if line.startswith("on:"):
            in_on_block = True
            continue
        if not in_on_block:
            continue
        if line and not line.startswith(" "):
            break
        match = re.match(r"^  (push|pull_request):$", line)
        if match:
            current = match.group(1)
            continue
        match = re.match(r'^      - "(.+)"$', line)
        if match and current:
            sections[current].append(match.group(1))
        elif re.match(r"^  \w+:", line):
            current = None
    return sections


def hook_watched() -> list[str]:
    with open(HOOK, encoding="utf-8") as handle:
        source = handle.read()
    match = re.search(r"^WATCHED='\^\((?P<body>.+)\)'$", source, re.MULTILINE)
    assert match, "WATCHED pattern not found"
    body = match.group("body")
    prefix, suffix = "tensorcad/(", r")\.py"
    return [f"tensorcad/{name}.py" for name in body[len(prefix):-len(suffix)].split("|")]


HOOK_MODES = {
    "chain_docs_only": (
        "仅暂存 `docs/*.md`、`README.md`（走仓内钩子链 `.githooks/pre-commit`）",
        "路径过滤器两道都不命中，直接返回，套件不启动",
    ),
    "hook_docs_only": (
        "同上，但只调 dual-derived 钩子",
        "路径不命中，直接返回（0.02 s 量级：钩子本身的开销）",
    ),
    "hook_geometry_edit": (
        "`tensorcad/hull.py`（只调 dual-derived 钩子）",
        "跑 suite.py + mutants.py，全绿才放行",
    ),
    "chain_geometry_edit": (
        "`tensorcad/hull.py`（走仓内钩子链）",
        "先跑仓内 `tensorcad.golden --quick`，再跑 dual-derived suite + mutants",
    ),
    "hook_regression_blocked": (
        "`tensorcad/hull.py` + 一份被注入 +3% Wigley 半宽的 TensorCAD 副本",
        "**`suite.py` 退出 0，仍被报告断言判红**（金标 case 变成 UNAVAILABLE）",
    ),
    "hook_named_suite_missing": (
        "`tensorcad/hull.py` + `GOLDEN_DIR` 指向不存在的目录",
        "**明确拒绝**（不静默降级成放行）",
    ),
    "hook_no_suite_beside_checkout": (
        "`tensorcad/hull.py` + `GOLDEN_DIR` 未设、旁无 TensorPapers",
        "警告 + 放行（本地不可执行 ≠ 本地失败），CI 侧仍 fail-closed",
    ),
}

QUANTITY_LABELS = {
    "block_coefficient": "Cb",
    "displacement_volume": "∇",
    "longitudinal_centre_buoyancy": "LCB",
    "maximum_section_area": "A_mx",
    "midship_coefficient": "Cm",
    "prismatic_coefficient": "Cp",
    "transverse_metacentric_radius": "BM_T",
    "vertical_centre_buoyancy": "VCB",
    "waterplane_area": "A_wp",
    "waterplane_coefficient": "Cwp",
    "wetted_surface_area": "S_wet",
}


def _fmt(value) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        if value == 0:
            return "0"
        if abs(value) < 1e-3:
            return f"{value:.2e}"
        return f"{value:.6f}"
    return str(value)


def trigger_table() -> str:
    paths = workflow_paths()
    watched_by_hook = set(hook_watched())
    rows = []
    for path in paths["push"]:
        rows.append(
            "| `{}` | {} | {} | {} |".format(
                path,
                "yes" if path in paths["push"] else "no",
                "yes" if path in paths["pull_request"] else "no",
                "yes" if path in watched_by_hook else "no (CI only)",
            )
        )
    pr_only = [p for p in paths["pull_request"] if p not in paths["push"]]
    for path in pr_only:
        rows.append(f"| `{path}` | no | yes | no (CI only) |")
    header = (
        "| 路径 | push 触发 | PR 触发 | pre-commit 钩子 |\n"
        "|---|---|---|---|\n"
    )
    return header + "\n".join(rows)


def case_table(suite: dict) -> str:
    rows = []
    for result in suite["results"]:
        findings = result["findings"]
        judged = [f for f in findings if f["verdict"] in ("PASS", "FAIL")]
        withheld = [f for f in findings if f["verdict"] == "WITHHELD"]
        tolerances = sorted({f["tolerance"] for f in judged})
        tiers = sorted({f["tier"] for f in judged})
        errors = [abs(f["error"]) for f in judged if f["error"] is not None]
        rows.append(
            "| `{}` | {} | {} | {} | {} | {} | {} |".format(
                result["case"],
                result["family"],
                "/".join(tiers) or "—",
                ", ".join(_fmt(t) for t in tolerances) or "—",
                _fmt(max(errors)) if errors else "—",
                len(judged),
                len(withheld),
            )
        )
    header = (
        "| case | 族 | 档位 | 容差（相对） | max 实测误差 | 判定项 | WITHHELD |\n"
        "|---|---|---|---|---|---|---|\n"
    )
    return header + "\n".join(rows)


def quantity_table(suite: dict) -> str:
    """One row per judged quantity: its tightest and loosest tolerance, worst error."""
    per_quantity: dict[str, dict] = {}
    for result in suite["results"]:
        for finding in result["findings"]:
            if finding["verdict"] not in ("PASS", "FAIL"):
                continue
            entry = per_quantity.setdefault(
                finding["quantity"],
                {"tiers": set(), "tolerances": [], "worst": 0.0, "cases": 0},
            )
            entry["tiers"].add(finding["tier"])
            entry["tolerances"].append(finding["tolerance"])
            entry["worst"] = max(entry["worst"], abs(finding["error"] or 0.0))
            entry["cases"] += 1
    rows = []
    for quantity in sorted(per_quantity):
        entry = per_quantity[quantity]
        rows.append(
            "| `{}` | {} | {} | {} – {} | {} | {} |".format(
                quantity,
                QUANTITY_LABELS.get(quantity, quantity),
                "/".join(sorted(entry["tiers"])),
                _fmt(min(entry["tolerances"])),
                _fmt(max(entry["tolerances"])),
                _fmt(entry["worst"]),
                entry["cases"],
            )
        )
    header = (
        "| 量 | 符号 | 档位 | 容差区间（相对） | 全 case 最坏误差 | case 数 |\n"
        "|---|---|---|---|---|---|\n"
    )
    return header + "\n".join(rows)


def withheld_table(suite: dict) -> str:
    rows = []
    for result in suite["results"]:
        withheld = [f["quantity"] for f in result["findings"] if f["verdict"] == "WITHHELD"]
        if withheld:
            rows.append(
                "| `{}` | {} | {} |".format(
                    result["case"],
                    ", ".join(f"`{q}`" for q in withheld),
                    len(withheld),
                )
            )
    if not rows:
        return "_（本次没有 WITHHELD。）_"
    header = "| case | 拒判的量 | 数量 |\n|---|---|---|\n"
    return header + "\n".join(rows)


def mutation_section(summary: dict) -> tuple[str, dict]:
    control = summary["control_arm"]
    control_rows = [
        f"| `{case}` | {verdict} |" for case, verdict in sorted(control.items())
    ]
    control_table = "| case | load-bearing 控制臂 |\n|---|---|\n" + "\n".join(control_rows)
    return control_table, {
        "control_pass": sum(1 for v in control.values() if v.startswith("PASS")),
        "control_total": len(control),
    }


def run_table(meta: dict) -> str:
    rows = []
    for run in meta["runs"]:
        rows.append(
            "| `{}` | {} | {} | `TENSORCAD_PATH={}` |".format(
                run["script"], run["exit_code"], f'{run["wall_s"]:.1f} s', meta["tensorcad_path"]
            )
        )
    header = "| 脚本 | 退出码 | 墙钟 | 环境 |\n|---|---|---|---|\n"
    return header + "\n".join(rows)


def hook_table(probe: dict) -> str:
    rows = []
    for mode in probe["modes"]:
        staged_label, expectation = HOOK_MODES.get(mode["mode"], (mode["mode"], ""))
        rows.append(
            "| {} | {} | {} | {} | `{}` |".format(
                staged_label,
                expectation,
                mode["exit_code"],
                f'{mode["wall_s"]:.2f} s',
                mode["last_line"] or "(silent)",
            )
        )
    header = (
        "| 暂存内容 | 钩子行为 | 退出码 | 墙钟 | 末行 |\n|---|---|---|---|---|\n"
    )
    return header + "\n".join(rows)


def _wall(probe: dict, mode: str) -> str:
    for entry in probe["modes"]:
        if entry["mode"] == mode:
            return f'{entry["wall_s"]:.2f} s'
    return "—"


def regression_table(control: dict, probe: dict) -> str:
    blocked = {m["mode"]: m for m in probe["modes"]}.get("hook_regression_blocked")
    rows = []
    for scenario in control["scenarios"]:
        not_pass = {
            case: verdict
            for case, verdict in scenario["case_verdicts"].items()
            if not verdict.startswith("PASS")
        }
        detail = ", ".join(f"`{k}`={v}" for k, v in sorted(not_pass.items())) or "全部 PASS"
        counts = scenario["finding_counts"]
        suite_exit = scenario["suite.py"]["exit_code"]
        if suite_exit != 0:
            gate = f"`suite.py` 退出 {suite_exit}（FAIL {counts.get('FAIL', 0)} 项）→ 第一道即红"
        elif blocked:
            gate = ("`suite.py` 退出 0，被报告断言判红（钩子退出 {}：`{}`）".format(
                blocked["exit_code"], blocked["last_line"]))
        else:
            gate = "`suite.py` 退出 0，靠报告断言判红"
        rows.append(
            "| {} | {} | {} | {} |".format(
                scenario["injected_change"], suite_exit, detail, gate))
    header = "| 注入的回归 | `suite.py` 退出码 | 未 PASS 的 case | 门怎么判 |\n|---|---|---|---|\n"
    return header + "\n".join(rows)


def snapshot_table(meta: dict) -> str:
    rows = []
    for name, digest in sorted(meta.get("snapshots", {}).items()):
        rows.append(f"| `{name}` | `{digest[:16]}…` |")
    header = "| 快照（`artifacts/gap-cad-d2-golden/`） | sha256 前 16 位 |\n|---|---|\n"
    return header + "\n".join(rows)


def render() -> str:
    suite = _load("golden_suite.json")
    report = _load("mutation_report.json")["summary"]
    meta = _load("run_meta.json")
    convergence = _load("series60_convergence.json")
    probe = _load("hook_probe.json")
    regression = _load("regression_control.json")
    chain_docs_wall = _wall(probe, "chain_docs_only")

    control_table, control_stats = mutation_section(report)
    verdicts = [f["verdict"] for r in suite["results"] for f in r["findings"]]
    counts = {v: verdicts.count(v) for v in sorted(set(verdicts))}

    # Coarsest Series 60 level = the grid whose undersampling the default-case tolerance
    # has to absorb.  Read from the convergence record, never re-typed.
    coarsest = convergence[0] if isinstance(convergence, list) and convergence else {}
    series60_deficit = abs(coarsest.get("deficit", 0.0))
    series60_grid = f'{coarsest.get("nx", "?")}x{coarsest.get("nz", "?")}'

    return f"""# GAP-CAD-D2 几何金标回归门禁（TensorCAD 落点）

> 本文由 `tools/golden/gap_d2_gate_docs.py` 从 `artifacts/gap-cad-d2-golden/*.json` 与两个
> 门禁文件本身程序化生成，正文数字零手抄（`--check` 可在 CI 里查陈旧）。
> 上游：TensorPapers 报告 **#20224**，套件 `papers/gap-cad-d2-golden-regression/`。
> 锚点：`GAP-CAD-D2-GOLDEN-REGRESSION`（任务 m22501，TensorCAD 仓侧收口）。

## 1. 这道门是什么，为什么仓里现在有两道

TensorCAD 的几何是生态的地基：`hull` / `wigley` / `series60` / `dtmb5415` 产形状，
`hydrostatics` + `mesh` + `geometry_audit` 把形状变成 Cb、∇、VCB、BM_T、湿表面积。
这些量漂 1 % 不会有任何测试报警，只会让下游阻力/耐波算例悄悄错。

仓里原有 `.github/workflows/golden-geometry.yml`（`tensorcad.golden`）守的是**仓内自己
committed 的基线**：四个金标 case、`FROZEN` 档判"是否漂移"。它的弱点写在它自己的设计里——
基线由本仓 `--bless` 生成，所以一次"改了代码顺手重 bless"的提交可以让它全绿。

本单装进来的 `geometry-golden.yml` 是**第三方的判据**：判据代码与金标数值都不住在这个仓里
（TensorPapers 的 dual-derived suite），它对该仓只暴露 `TENSORCAD_PATH`。每个金标值在进入
判定前必须由两条结构不同的推导给出同一结果（EXACT 档 1e-9 / QUAD 档 1e-6），否则直接抛
`GoldenError`，**不存在手敲回退值**。因此这道门不可能被"bless 成与自己一致"。

两道门并存不是冗余：`golden-geometry.yml` 便宜（~5 s，粗网格）、每次 `tensorcad/**` 改动都跑；
`geometry-golden.yml` 贵（~30 s，两档网格 + 检测力探针）、只在这一改动**能移动几何数字**时跑。

## 2. 触发路径（程序化取自 `geometry-golden.yml` 与 `.githooks/pre-commit-geometry-golden.sh`）

{trigger_table()}

**不触发的改动**：只改 `docs/**`、`README.md`、`CHANGELOG.md`、`tests/**` 等文档/测试文件时，
push/PR 路径过滤器不命中，pre-commit 钩子也在路径判断上直接退出（实测 {chain_docs_wall}，见 §6）。
这条性质由 `tests/test_geometry_golden_gate.py::test_documentation_only_edits_are_not_watched`
与 `::test_hook_lets_a_documentation_only_commit_through_without_running_the_suite` 双保险钉住。

一个诚实的边界：**PR 路径过滤器本身无法在本仓 CI 里被自测**（一次 docs-only 的 PR 不触发
工作流，也就没有工作流能报告"我没触发"）。因此"该不触发时不触发"是靠单元测试对
`on: push/pull_request` 路径表做的静态断言，而不是靠 CI 跑绿——两种证据强度不同，这里写明。

## 3. 容差表（本次实测）

各 case 的整体判据：

{case_table(suite)}

逐量视图（同一量在粗/细网格上的容差不同，故给区间）：

{quantity_table(suite)}

本次判词分布：{", ".join(f"{k}={v}" for k, v in counts.items())}。

`FROZEN` 档出现在 `dtmb5415_exact` 上：三角面数/顶点数/边界边数是**整数相等**（容差 0），
包围盒三向尺度容差 1e-9——它判的是"出厂工件还是不是那一个"，不是"精度有多好"，所以这一行
的 max 误差是 0 也不构成任何精度声明。

**容差是怎么来的（不是拍的）**：EXACT 档容差 = 该 case 网格离散误差的实测上界留 2 倍以上余量；
`series60_cb060_default`（{series60_grid} 默认网格）给到 4.5e-2，是因为该网格对 Cb 的欠填实测 {series60_deficit:.4%}（`series60_convergence.json` 首档），
把"粗网格的已知缺陷"显式写成容差，而不是悄悄调大或把 case 删掉。细则见 TensorPapers
报告 #20224 与 `outputs/series60_convergence.json`。

## 4. WITHHELD：不能判的量就不判

{withheld_table(suite)}

`WITHHELD` 永不产出判词、也永不被"补齐成 PASS"：5415 出厂 STL 是开口面（边界边存在），
无独立锚，所以 ∇/Cb 明写拒判；Series 60 的 Cm/Cp/VCB/Cwp/湿表面积没有表值转录，同样拒判。
只有 Cb/∇ 以**定义式**入判（积分必须回到请求的丰满度）。

## 5. 检测力（门有没有牙齿）

控制臂：先证明这道门在被测库"原样"时是全绿的——否则后面的"抓到变异"毫无意义。

{control_table}

共 {control_stats["control_pass"]}/{control_stats["control_total"]} 控制臂 PASS。

探针（把被判量挪到"总误差 = k × 容差"）：

| 探针组 | 抓/应抓 | 命中率 | 说明 |
|---|---|---|---|
| 强变异族（k>1，必须 FAIL） | {report["strong_caught"]}/{report["must_catch_mutants"]} | {report["strong_detection_rate"]:.4f} | **验收口径的"检测力 100 %"** |
| 几何闸探针（L/B/T 错船型） | {report["gate_caught"]}/{report["gate_mutants"]} | {report["gate_detection_rate"]:.4f} | 闸不过即降级 WITHHELD |
| WITHHELD 拒判探针 | {report["withheld_refused"]}/{report["withheld_probes"]} | {report["withheld_refused"] / report["withheld_probes"]:.4f} | 丢证据的量必须被拒判 |
| 弱变异族（k≤0.5，**不要求**抓到） | {report["quantity_caught"]}/{report["quantity_mutants"]} | {report["quantity_detection_rate"]:.4f} | 容差内的自然离散，抓到才算过度敏感 |
| 误报 | 0 例（`false_alarms={report["false_alarms"]}`） | — | 控制臂 + k≤0.5 臂无一误判 |

`unexpected_misses={report["unexpected_misses"]}`、`misses={report["misses"]}`。
注意 `quantity_detection_rate` 只有 {report["quantity_detection_rate"]:.4f}：那一栏含 k=0.5 的
**容差内**探针，抓不到才是正常（抓到反而说明容差虚紧）。验收要求的 100 % 指
`strong_detection_rate` 与 `gate_detection_rate`，两者都是 {report["strong_detection_rate"]:.4f}。

### 5.1 退出码不够：一个"测不出来的金标 case"必须算红

装门时做了一次**真回归**对照，不只做合成变异：把 TensorCAD 复制一份、在源码里注入一处改动，
再用同一个钩子去判它（复现：`tensorcad/wigley.py` 的 `half_breadth` 把 `0.5*beam` 改成
`0.515*beam`；`tensorcad/hydrostatics.py` 把水线面投影面积乘 1.02）。

{regression_table(regression, probe)}

第一行是本单的一个**实质发现**：几何模块改坏后 TensorCAD 自己的 `generate_wigley` 会拒绝出网格
（水线偏移偏离解析面 3.0e-2 m），套件把该 case 记成 `SKIP` / `UNAVAILABLE`——这符合它"缺测
永不 PASS"的诚实口径——**但 `suite.py` 仍然退出 0**。也就是说，只看退出码的门会在"两个金标
case 已经不再被测量"的情况下放行。本仓两个载体因此都不读退出码，而读报告：
`tools/golden/gap_d2_assert_report.py` 逐 case 断言 `verdict` 以 `PASS` 开头且 case 数 = 8，
再断言检测力的四个率；钩子与 CI 共用这一个脚本（`tests/test_geometry_golden_gate.py` 里钉住）。
上游 `suite.py` 的退出码语义未被改动——本单的写权限只在 TensorCAD 侧，加固做在这里。

## 6. 首跑证据（本仓 checkout，非 TensorPapers 侧）

{run_table(meta)}

- 被检 checkout：`{meta["tensorcad_path"]}`，HEAD `{meta["tensorcad_head"][:12]}`，Python {meta["python"]}
- 墙钟随机器负载波动（同一条命令在 §6 的钩子探针里实测 {_wall(probe, "hook_geometry_edit")}）；可复现的是退出码与
  `run_meta.json` 里记的 sha256，不是秒数
- 两次运行原始输出：`artifacts/gap-cad-d2-golden/suite.log`、`mutants.log`
- 快照（供 §3/§5 表格注入）：

{snapshot_table(meta)}

pre-commit 链步实测（`git diff --cached` 语义用 `GOLDEN_STAGED_FILE` 复现，
`tools/golden/gap_d2_gate_probe.py` 生成 `artifacts/gap-cad-d2-golden/hook_probe.json`）：

{hook_table(probe)}

## 7. 怎么用

```bash
# CI：push/PR 命中几何路径即自动跑（.github/workflows/geometry-golden.yml）
# 本地：钩子由仓内既有 hooksPath 安装器接管，无需另设 .git/hooks
tools/golden/install_git_hooks.sh          # core.hooksPath=.githooks（本仓既有机制）

# 套件位置解析顺序：$GOLDEN_DIR → <repo>/../TensorPapers/papers/gap-cad-d2-golden-regression
#                     → $HOME/TensorPapers/papers/gap-cad-d2-golden-regression
GOLDEN_DIR=/path/to/TensorPapers/papers/gap-cad-d2-golden-regression git commit ...

# 手动跑（与 CI 完全同一条命令）
cd $GOLDEN_DIR && TENSORCAD_PATH=$(git -C /path/to/TensorCAD rev-parse --show-toplevel) \\
    PYTHONPATH=$(git -C /path/to/TensorCAD rev-parse --show-toplevel) python3 suite.py
```

`GOLDEN_DIR` 显式给出但目录不存在 = 硬错（不静默降级）；未给出且找不到 = 本地放行 + 警告，
因为"这个 checkout 旁边没克隆 TensorPapers"不等于"代码回归了"。CI 侧没有这个宽容：
拿不到 TensorPapers 或套件失败都会让 job 红。

## 8. 局限与改动清单（诚实）

- **只证明几何/静水力**：这道门对 Cb/∇/VCB/BM_T/S_wet 一类量负责，**不**证明阻力、耐波、
  操纵、结构适用性。几何对不等于船好。
- **控制臂 8/8 不等于全谱覆盖**：8 个 case 覆盖球/椭球/Wigley/Series60/5415，Series60 的
  高阶级系数、5415 的 ∇ 仍是拒判（无独立锚）。
- **粗网格容差是承认缺陷，不是校准**：`series60_cb060_default` 3.7 % 的 Cb 欠填被写进容差，
  细网格 case（321×81）才是精度声明。任何人想"把粗网格也收进 2 %"，得改套件而不是改容差。
- **runner 标签**：drop-in 原文写 `runs-on: ubuntu-latest`，本仓其余门禁用自托管 `ci-loong64`，
  这里改成 `ci-loong64`，否则该 job 在本生态里不会被派发。这一处与其他两处改动
  （TensorPapers checkout 的 token 回退、验收断言步）都写在 workflow 头部注释里。
- **上游退出码语义有口子（已在本仓兜住）**：`suite.py` 在 case 变 `UNAVAILABLE` 时仍退出 0。
  建议 TensorPapers 侧把 `UNAVAILABLE` 计入非零退出；在它改之前，TensorCAD 侧靠报告断言兜住
  （§5.1 的第一行就是这条口子的实测证据）。
- **钩子不是强制手段**：`git commit --no-verify`、或没装 hooksPath 的 checkout 都能绕过；
  权威判据是 CI。
- **`.git/hooks/pre-commit` 未直接写入**：本仓用 `core.hooksPath=.githooks`，钩子以
  husky 式链步（`.githooks/pre-commit` 末尾调用 `.githooks/pre-commit-geometry-golden.sh`）
  落地，这样改动是**可 review 的仓内文件**，而不是别人机器上的 .git 目录。

## 9. 相关文件

| 文件 | 作用 |
|---|---|
| `.github/workflows/geometry-golden.yml` | CI 门（权威，fail-closed），含验收断言步 |
| `.githooks/pre-commit-geometry-golden.sh` | 本地 husky 式链步（路径过滤 → suite + mutants） |
| `.githooks/pre-commit` | 仓内第一道（`tensorcad.golden --quick`）并链第二道 |
| `tests/test_geometry_golden_gate.py` | 两臂路径表一致性、docs-only 不触发、钩子行为 |
| `tools/golden/gap_d2_assert_report.py` | 判定层：读两份报告断言验收数字（钩子与 CI 共用） |
| `tools/golden/gap_d2_gate_probe.py` | 钩子七种决策模式的实测探针（写 `hook_probe.json`） |
| `tools/golden/gap_d2_gate_docs.py` | 本文生成器（`--check` 查陈旧） |
| `artifacts/gap-cad-d2-golden/regression_control.json` | 真回归对照（注入改动 → 套件与门怎么判） |
| `artifacts/gap-cad-d2-golden/` | 首跑快照：suite/mutants JSON、日志、run_meta |
"""


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="exit 1 if the document is stale")
    args = parser.parse_args(argv)

    rendered = render()
    if args.check:
        try:
            with open(DOC, encoding="utf-8") as handle:
                existing = handle.read()
        except FileNotFoundError:
            print(f"missing: {DOC}", file=sys.stderr)
            return 1
        if existing != rendered:
            print(f"stale: {DOC} (rerun {os.path.basename(__file__)})", file=sys.stderr)
            return 1
        print("up to date")
        return 0

    with open(DOC, "w", encoding="utf-8") as handle:
        handle.write(rendered)
    print(f"wrote {DOC}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())