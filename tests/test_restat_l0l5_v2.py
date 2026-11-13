"""AUDIT-RESTAT v2 口径逻辑单测（脱库：只测纯函数，不连 PG）。

对应实现：tools/audit_restat/restat_v2.py
被测对象是「口径」本身——锚点抽取、scope 解析、claim-gate 严格性、定义敏感带、
溯源漏斗。这些函数一旦漂移，报告里的所有数字都会跟着漂，故用脱库单测钉住。
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "audit_restat"))

import restat_v2 as R  # noqa: E402


# ---------------------------------------------------------------- 锚点抽取

@pytest.mark.parametrize("goal,expected", [
    ("【锚点=AUDIT-RESTAT】【自找活·tp-glm51-2】正文", "AUDIT-RESTAT"),
    ("【锚点=CASE-LBM-SQUARE-CYL-NORBERG】\n多行正文", "CASE-LBM-SQUARE-CYL-NORBERG"),
    ("前言\n【锚点=GAP-ICE-CRUSH-ISO】后有文字", "GAP-ICE-CRUSH-ISO"),
    ("无锚点的历史单", None),
    (None, None),
    ("", None),
])
def test_extract_anchor(goal, expected):
    assert R.extract_anchor(goal) == expected


@pytest.mark.parametrize("goal", [
    "【锚点=SLUG】模板未替换",
    "【锚点=NEW:】空理由",
    "【锚点=NEW:理由】占位理由",
])
def test_extract_anchor_rejects_placeholders(goal):
    """占位串必须抽成 None：否则 770 张模板单会被当成真实能力单元，锚点数虚高。"""
    assert R.extract_anchor(goal) is None


# ---------------------------------------------------------------- scope 解析

def test_scope_anchor_parses_canonical_form():
    assert R.scope_anchor("m17013 KPI-LBM-VERIFIED-CAPS") == "KPI-LBM-VERIFIED-CAPS"
    assert R.scope_anchor("m16243 CASE-DEM-SHEAR-CAMPBELL-ELASTIC-LIMIT") == \
        "CASE-DEM-SHEAR-CAMPBELL-ELASTIC-LIMIT"


@pytest.mark.parametrize("scope", [
    "hydrostatics_analytic",
    "docs/ontology v1 + card schemas + validate_ontology_v1 (m8815)",
    "fluctuating_hydro (PR#433) @HEAD-red",
    "",
    None,
])
def test_scope_anchor_rejects_freeform(scope):
    assert R.scope_anchor(scope) is None


# ---------------------------------------------------------------- claim vs gate

@pytest.mark.parametrize("claim,gate,expected", [
    (4, 4, "match"),
    (4, 1, "overclaim"),
    (2, 3, "underclaim"),
    (3, None, "ungated"),
    (1, 0, "overclaim"),
    (0, 0, "match"),
])
def test_strictness(claim, gate, expected):
    assert R.strictness(claim, gate) == expected


def test_crosstab_counts_and_ungated_column():
    rows = [
        {"level": 2, "gate_level": 2},
        {"level": 2, "gate_level": 1},
        {"level": 3, "gate_level": None},
        {"level": 0, "gate_level": 0},
        {"level": 9, "gate_level": 9},   # 越界 claim 被丢弃
    ]
    ct = R.crosstab(rows)
    assert ct["claim L2"]["2"] == 1
    assert ct["claim L2"]["1"] == 1
    assert ct["claim L3"]["ungated"] == 1
    assert ct["claim L0"]["0"] == 1
    assert set(ct) == {f"claim L{i}" for i in range(6)}


# ---------------------------------------------------------------- 比率与带

def test_ratio_pct_handles_zero_denominator():
    assert R.ratio(3, 0) is None
    assert R.pct(1, 0) is None
    assert R.pct(1, 4) == 25.0
    assert R.pct2(1, 3000) == 0.03


def test_band_is_min_max_not_mean():
    b = R.band([2313, 2718, 2856, None])
    assert b == {"min": 2313, "max": 2856, "n_definitions": 3}
    assert R.band([]) is None


# ---------------------------------------------------------------- Tier-B 定义带

def _m(files, pr_state="merged"):
    return {"changed_files": files, "pr_state": pr_state}


def test_tier_b_definitions_and_band():
    missions = [
        _m(["tensorcad/hull.py"]),                       # A,B,C 计数
        _m(["tests/test_hull.py"]),                       # 仅 A,C（测试不算实现文件）
        _m(["docs/hull.md"]),                             # 仅 A,C
        _m([], "merged"),                                 # A 计（merged 但无文件），C 不计
        _m(["tensorcad/solver.cpp"], "open"),             # 仅 D（未入库实现）
        _m(["tensorcad/verify_report.py"], "merged"),     # 仅 A,C（helper 命名不算实现）
    ]
    out = R.tier_b_strict_bands(missions)
    d = out["definitions"]
    assert d["A_merged_any_files"] == 5      # 6 条里 5 条 merged（solver.cpp 那条是 open）
    assert d["B_merged_impl_file"] == 1      # 仅 tensorcad/hull.py 算实现文件
    assert d["C_merged_nonempty_changed_files"] == 4
    assert d["D_impl_file_any_pr_state"] == 2  # 不看合并：hull.py + solver.cpp
    # 带只取 A/B/C（要求合并），不含宽松上界 D
    assert out["band"] == {"min": 1, "max": 5, "n_definitions": 3}
    assert out["loose_upper_bound"] == 2


# ---------------------------------------------------------------- 溯源漏斗

def test_traceability_funnel_counts():
    ledger = [
        {"slug": "A-ONE", "status": "verified"},
        {"slug": "A-TWO", "status": "done"},
        {"slug": "A-THREE", "status": "open"},
    ]
    verdicts = [
        {"repo": "X", "scope": "m1 A-ONE", "level": 4, "gate_level": 4, "verdict_by": "p"},
        {"repo": "X", "scope": "m2 A-TWO", "level": 2, "gate_level": 1, "verdict_by": "p"},
        {"repo": "X", "scope": "freeform_scope", "level": 1, "gate_level": 1, "verdict_by": "p"},
    ]
    missions = [
        {"goal": "【锚点=A-ONE】", "pr_state": "merged", "changed_files": ["tensorcad/a.py"]},
        {"goal": "【锚点=A-ONE】", "pr_state": "open", "changed_files": []},
        {"goal": "【锚点=SLUG】", "pr_state": "merged", "changed_files": []},
        {"goal": "无锚点", "pr_state": None, "changed_files": []},
    ]
    tr = R.traceability(ledger, verdicts, missions)
    assert tr["missions_succeeded"] == 4
    assert tr["missions_merged"] == 2
    assert tr["missions_with_anchor"] == 2
    assert tr["missions_anchor_placeholder"] == 1
    assert tr["distinct_anchors_in_missions"] == 1
    assert tr["tasks_per_anchor"] == 2.0
    assert tr["verdict_scopes_parsable"] == 2
    assert tr["verdict_scopes_freeform"] == 1
    assert tr["verdict_scopes_mapped_to_ledger"] == 2
    assert tr["verdict_scopes_unmapped_to_ledger"] == 0
    assert tr["succeeded_to_ledger_pct"] == 75.0
    assert tr["ledger_verified"] == 1


def test_traceability_flags_unmapped_scope_anchor():
    ledger = [{"slug": "A-ONE", "status": "done"}]
    verdicts = [{"repo": "X", "scope": "m7 GHOST-ANCHOR", "level": 2, "gate_level": 2,
                 "verdict_by": "p"}]
    tr = R.traceability(ledger, verdicts, [])
    assert tr["verdict_scopes_unmapped_to_ledger"] == 1
    assert tr["unmapped_sample"] == ["GHOST-ANCHOR"]


def test_traceability_survives_empty_inputs():
    """库空时不得崩：报告工具要能在半初始化环境跑出 0 而不是 traceback。"""
    tr = R.traceability([], [], [])
    assert tr["missions_succeeded"] == 0
    assert tr["tasks_per_anchor"] is None
    assert tr["succeeded_to_ledger_pct"] is None


# ---------------------------------------------------------------- Tier-A 汇总

def test_tier_a_stats_by_producer_split():
    verdicts = [
        {"repo": "R", "scope": "m1 A", "level": 4, "gate_level": 1, "verdict_by": "strict"},
        {"repo": "R", "scope": "m2 B", "level": 3, "gate_level": 3, "verdict_by": "lenient"},
        {"repo": "R", "scope": "m3 C", "level": 2, "gate_level": None, "verdict_by": "lenient"},
    ]
    ta = R.tier_a_stats(verdicts)
    assert ta["n_verdicts"] == 3
    assert ta["n_scopes"] == 3
    assert ta["gate2"] == 1
    assert ta["overclaim"] == 1
    assert ta["by_producer"]["strict"]["overclaim_rate_pct"] == 100.0
    assert ta["by_producer"]["lenient"]["gate2_rate_pct"] == 50.0
    assert ta["by_producer"]["lenient"]["ungated"] == 1