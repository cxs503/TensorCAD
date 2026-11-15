"""Evidence bundle for the large-geometry performance benchmark.

``docs/perf_scale_10m.md`` reports timings for 10^6--10^7 triangle geometry.
Those numbers used to be transcribed by hand out of ``artifacts/perf_scale/``,
and ``artifacts/`` is stripped by the merge gate -- so the published table was
neither reproducible from the repository nor checkable against anything.  This
module closes that hole: the result tables in the document are *rendered* from
the committed measurement JSON, and :func:`check_document` fails when the text
drifts from the artifact.

The contract is deliberately dumb and stdlib-only (TensorCAD is zero-dependency):

* a table lives between ``<!-- perf-evidence:<name> begin -->`` and
  ``<!-- perf-evidence:<name> end -->`` markers;
* the block body is exactly ``render_block(name, artifacts)``;
* :func:`check_document` compares the two strings, so a stale or hand-edited
  number is a test failure rather than a silent claim.

Usage::

    python scripts/perf_scale_doc.py --write   # re-render from artifacts
    python scripts/perf_scale_doc.py --check   # CI: doc matches artifacts
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, Sequence

#: Measured sub-tasks, in display order, with the axis they scale in.
MESH_ORDER = ("load", "point_locate", "point_locate_accel", "volume")
#: Sub-tasks on the multi-body (fleet) 3-point ladder, in display order.
FLEET3PT_ORDER = ("load", "fleet_assemble", "point_locate_accel", "volume")
AXIS_LABEL = {
    "load": ("三角面", "tri/s"),
    "point_locate": ("三角面", "query/s"),
    "point_locate_accel": ("三角面", "query/s"),
    "volume": ("三角面", "tri/s"),
    "voxelize": ("体素", "cell/s"),
    "fleet_assemble": ("三角面", "tri/s"),
}

DASH = "—"

#: Wall-clock budget used to turn a measured rate into a "how far can one
#: machine go" number.  The rule is stated in the document so the cluster
#: boundary in acceptance item 3 is reproducible instead of asserted.
BUDGET_SECONDS = 60.0


def marker_begin(name: str) -> str:
    return f"<!-- perf-evidence:{name} begin -->"


def marker_end(name: str) -> str:
    return f"<!-- perf-evidence:{name} end -->"


def _num(value: Any, spec: str) -> str:
    """Format a measurement, rendering ``None``/missing as an em dash."""
    if value is None:
        return DASH
    try:
        return format(float(value), spec)
    except (TypeError, ValueError):  # pragma: no cover - defensive
        return DASH


def _int(value: Any) -> str:
    if value is None:
        return DASH
    return f"{int(value):,}"


def _results(artifact: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = list(artifact.get("results") or [])
    if not rows:
        raise ValueError("artifact has no results")
    return rows


def _measured(artifact: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [r for r in _results(artifact) if r.get("status") == "measured"]


def render_mesh_table(artifact: Mapping[str, Any]) -> str:
    """Markdown table of the mesh sub-tasks, one row per measurement.

    Columns answer the acceptance questions directly: how long, how fast, how
    much did the acceleration structure cost to build, and how much memory did
    the process hold.  ``peak_rss_mb`` is a process-wide high-water mark
    (``ru_maxrss``), so it is labelled as such rather than read per rung.
    """
    rows = _results(artifact)
    order = {name: i for i, name in enumerate(MESH_ORDER)}
    rows = sorted(
        rows,
        key=lambda r: (order.get(r.get("subtask"), len(order)), r.get("triangles") or 0),
    )
    lines = [
        "| 子任务 | 三角面 | 秒 | 吞吐 | 建索引 s | 索引 MB | 峰值 RSS MiB |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        subtask = r.get("subtask", "?")
        if r.get("copies", 1) > 1:
            subtask = f"{subtask} ×{r['copies']}"
        _, unit = AXIS_LABEL.get(r.get("subtask"), ("", ""))
        if r.get("status") != "measured":
            cells = [
                f"**{subtask}**",
                _int(r.get("triangles")),
                "not_run",
                r.get("note") or "",
                DASH,
                DASH,
                _num(r.get("peak_rss_mb"), ",.0f"),
            ]
        else:
            index_mb = (r.get("index_bytes") or 0) / (1024.0 * 1024.0)
            cells = [
                f"**{subtask}**" if r.get("subtask") == "point_locate_accel" else subtask,
                _int(r.get("triangles")),
                _num(r.get("seconds"), ".3f"),
                f"{_num(r.get('throughput'), ',.0f')} {unit}".strip(),
                _num(r.get("build_seconds"), ".2f") if r.get("build_seconds") else DASH,
                _num(index_mb, ".1f") if index_mb else DASH,
                _num(r.get("peak_rss_mb"), ",.0f"),
            ]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def render_voxel_table(artifact: Mapping[str, Any]) -> str:
    """Markdown table of the voxelisation ladder (cells, seconds, cells/s, RSS)."""
    rows = sorted(_results(artifact), key=lambda r: r.get("grid_cells") or 0)
    lines = [
        "| 网格 | 体素 | 秒 | 吞吐 cell/s | 峰值 RSS MiB |",
        "|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        measured = r.get("status") == "measured"
        lines.append(
            "| "
            + " | ".join([
                f"{r.get('nx')}³",
                _int(r.get("grid_cells")),
                _num(r.get("seconds"), ".3f") if measured else "not_run",
                _num(r.get("throughput"), ",.0f") if measured else (r.get("note") or ""),
                _num(r.get("peak_rss_mb"), ",.0f"),
            ])
            + " |"
        )
    return "\n".join(lines)


def render_budget_table(
    mesh: Mapping[str, Any],
    voxel: Mapping[str, Any],
    budget_seconds: float = BUDGET_SECONDS,
) -> str:
    """One row per (sub-task, body count): the largest measured rung and the
    size the *measured* rate would reach in ``budget_seconds``.

    This is the reproducible half of "at what scale does this need a cluster":
    the extrapolation is multiplication of a measured rate by a stated budget,
    so anyone can re-derive it from the artifact (no hand-tuned thresholds).
    """
    groups: dict[tuple[str, int], list[dict[str, Any]]] = {}
    for artifact in (mesh, voxel):
        for r in _measured(artifact):
            if not r.get("seconds") or not r.get("throughput"):
                continue
            groups.setdefault((r["subtask"], r.get("copies", 1)), []).append(r)

    lines = [
        f"| 子任务 | 规模轴 | 实测最大档 | 秒 | 实测速率 | {budget_seconds:.0f} s 预算等价规模 | 峰值 RSS MiB |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    order = {name: i for i, name in enumerate(MESH_ORDER + ("fleet_assemble", "voxelize"))}

    def biggest(rows: list[dict[str, Any]]) -> dict[str, Any]:
        return max(rows, key=lambda r: (r.get("triangles") or 0, r.get("grid_cells") or 0))

    for key in sorted(groups, key=lambda k: (order.get(k[0], len(order)), k[1])):
        subtask, copies = key
        r = biggest(groups[key])
        axis_name, unit = AXIS_LABEL.get(subtask, ("规模", ""))
        equivalent = r["throughput"] * budget_seconds
        label = subtask if copies == 1 else f"{subtask} ×{copies}"
        size = r.get("triangles") or r.get("grid_cells") or 0
        lines.append(
            "| "
            + " | ".join([
                label,
                axis_name,
                _int(size),
                _num(r.get("seconds"), ".3f"),
                f"{_num(r.get('throughput'), ',.0f')} {unit}".strip(),
                _int(equivalent),
                _num(r.get("peak_rss_mb"), ",.0f"),
            ])
            + " |"
        )
    return "\n".join(lines)


def render_fleet3pt_table(artifact: Mapping[str, Any]) -> str:
    """Markdown table of the multi-body (fleet) 3-point ladder.

    One row per (sub-task, body count, size).  The 10^6--10^7 band is only
    reachable with disjoint copies of one generated body, so the body count is
    part of the row identity: ``load`` is always a single body (it *is* the
    generator), ``fleet_assemble``/``point_locate_accel``/``volume`` run on 2
    and 3 bodies.  ``peak_rss_mb`` is the process high-water mark (labelled as
    such), ``rss_mb``/``mesh_mb`` are the per-rung live numbers.
    """
    rows = _results(artifact)
    order = {name: i for i, name in enumerate(FLEET3PT_ORDER)}
    rows = sorted(
        rows,
        key=lambda r: (order.get(r.get("subtask"), len(order)),
                       r.get("copies") or 1, r.get("triangles") or 0),
    )
    lines = [
        "| 子任务 | 体数 | 三角面 | 秒 | 吞吐 | 建索引 s | 索引 MiB | 活 RSS MiB | 峰值 RSS MiB |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        subtask = r.get("subtask", "?")
        _, unit = AXIS_LABEL.get(subtask, ("", ""))
        measured = r.get("status") == "measured"
        index_mb = (r.get("index_bytes") or 0) / (1024.0 * 1024.0)
        lines.append(
            "| "
            + " | ".join([
                f"**{subtask}**" if subtask == "point_locate_accel" else subtask,
                str(r.get("copies") or 1),
                _int(r.get("triangles")),
                _num(r.get("seconds"), ".3f") if measured else "not_run",
                f"{_num(r.get('throughput'), ',.0f')} {unit}".strip() if measured
                else (r.get("note") or ""),
                _num(r.get("build_seconds"), ".2f") if r.get("build_seconds") else DASH,
                _num(index_mb, ".1f") if index_mb else DASH,
                _num(r.get("rss_mb"), ",.0f"),
                _num(r.get("peak_rss_mb"), ",.0f"),
            ])
            + " |"
        )
    return "\n".join(lines)


def render_fleet3pt_gates(artifact: Mapping[str, Any]) -> str:
    """Markdown table of the complexity (exponent) gates on the fleet ladder.

    A gate that could not be evaluated is printed as ``not_run`` -- never as a
    pass -- and the footer repeats the three-state ``verdict`` so a reader
    cannot mistake "no data" for "within limit".  Exponent and R^2 are parsed
    out of the gate ``detail`` string with a regex instead of a human reading
    them, and a gate whose detail does not parse is rendered as ``—`` rather
    than silently mis-filled.
    """
    import re

    lines = [
        "| 门 | 判定 | 拟合指数 | 限值 | R² | 规模点数 | 细节 |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for g in artifact.get("gates", []):
        status = g.get("status")
        detail = g.get("detail") or ""
        exponent = re.search(r"exponent=([-0-9.]+)", detail)
        limit = re.search(r"limit ([0-9.]+)", detail)
        r2 = re.search(r"R\^2=([0-9.]+)", detail)
        count = re.search(r"(\d+) sizes", detail)
        if status == "not_run":
            mark = "not_run（未检查）"
        else:
            mark = "pass" if g.get("passed") else "**RED**"
        lines.append(
            "| "
            + " | ".join([
                f"`{g.get('name')}`",
                mark,
                exponent.group(1) if exponent else DASH,
                limit.group(1) if limit else DASH,
                r2.group(1) if r2 else DASH,
                count.group(1) if count else DASH,
                detail,
            ])
            + " |"
        )
    checks = artifact.get("volume_checks") or []
    footer = (
        f"\n\nverdict=**{artifact.get('verdict')}**，gates_checked={artifact.get('gates_checked')}，"
        f"gates_not_run={artifact.get('gates_not_run')}"
    )
    if checks:
        worst = max(checks, key=lambda c: abs(c.get("rel_err") or 0.0))
        footer += (
            f"；体积校验 {sum(1 for c in checks if c.get('ok'))}/{len(checks)} 通过"
            f"（最差 {worst.get('rel_err'):+.4%} @ {_int(worst.get('triangles'))} 三角面，"
            f"判据 3%）"
        )
    return "\n".join(lines) + footer


def render_baseline_gate_table(baseline: Mapping[str, Any]) -> str:
    """Markdown table of the committed quick-ladder wall-clock baseline.

    This is the reference every 2x regression gate compares against, so it is
    rendered from ``benchmarks/perf_scale/baseline.json`` rather than copied
    into prose: the gate keys are built with
    :func:`tensorcad.perf_benchmark.measurement_key`, i.e. exactly the keys the
    gate lookup uses, and the RED threshold column uses the module default
    factor instead of a hand-typed ``2``.
    """
    from .perf_benchmark import DEFAULT_GATE_FACTOR, measurement_key

    rows = _results(baseline)
    order = {name: i for i, name in enumerate(MESH_ORDER + ("voxelize",))}
    rows = sorted(
        rows,
        key=lambda r: (order.get(r.get("subtask"), len(order)),
                       r.get("nx") or 0, r.get("nz") or 0),
    )
    factor = DEFAULT_GATE_FACTOR
    lines = [
        f"| 门键 | 规模 | 基线秒 | RED 阈值 {factor:g}× | 基线吞吐 |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in rows:
        subtask = r.get("subtask", "?")
        axis_name, unit = AXIS_LABEL.get(subtask, ("规模", ""))
        size = r.get("triangles") or r.get("grid_cells") or 0
        seconds = r.get("seconds")
        key = measurement_key(subtask, int(r.get("nx") or 0), int(r.get("nz") or 0),
                              int(r.get("copies") or 1))
        lines.append(
            "| "
            + " | ".join([
                f"`{key}`",
                f"{axis_name} {_int(size)}",
                _num(seconds, ".4f"),
                _num(None if seconds is None else seconds * factor, ".4f"),
                f"{_num(r.get('throughput'), ',.0f')} {unit}".strip(),
            ])
            + " |"
        )
    env = baseline.get("environment") or {}
    lines.append("")
    lines.append(
        f"基线记录 {len(rows)} 条；环境 `machine={env.get('machine')}`、"
        f"`cpu_count={env.get('cpu_count')}`、`python={env.get('python')}`；"
        f"口径：{baseline.get('baseline_policy')}"
    )
    return "\n".join(lines)


def render_code_facts() -> str:
    """Markdown table of the benchmark's tuning constants, read from the code.

    A document that quotes a constant is a second source of truth; this table
    removes the second source.  Values are read from
    :mod:`tensorcad.perf_benchmark` at render time, so changing a threshold
    without re-rendering is a test failure.
    """
    from . import perf_benchmark as pb

    def ladder(pairs: tuple[tuple[int, int], ...]) -> str:
        return " → ".join(f"{nx}×{nz}" for nx, nz in pairs)

    facts: list[tuple[str, str, str]] = [
        ("`CLUSTER_THRESHOLD`", f"{pb.CLUSTER_THRESHOLD:,}",
         "朴素路径的单次规模告警线（三角面）；不是实测集群边界"),
        ("`DEFAULT_GATE_FACTOR`", f"{pb.DEFAULT_GATE_FACTOR:g}",
         "墙钟回归门 RED 阈值（`--regression-factor` 的默认值）"),
        ("`DEFAULT_MAX_NAIVE_TRIANGLES`", f"{pb.DEFAULT_MAX_NAIVE_TRIANGLES:,}",
         "超过即跳过 O(N·Q) 暴力参考，标 `not_run`"),
        ("`MIN_SCALING_POINTS`", f"{pb.MIN_SCALING_POINTS}",
         "规模点数不足则指数门 `not_run`（不是 pass）"),
        ("`MESH_CACHE_SIZE`", f"{pb.MESH_CACHE_SIZE}", "阶梯每档只生成一次网格的缓存上限"),
        ("`DEFAULT_LADDER`", ladder(pb.DEFAULT_LADDER),
         "full 阶梯（顶端受 `HullParameters` 上限约束）"),
        ("`QUICK_LADDER`", ladder(pb.QUICK_LADDER), "quick 阶梯，即 CI 基线覆盖的范围"),
        ("`VOXEL_LADDER`", " → ".join(f"{g}³" for g in pb.VOXEL_LADDER), "体素阶梯"),
        ("`SCALING_LIMITS`",
         "; ".join(f"`{k}`≤{v:g}" for k, v in sorted(pb.SCALING_LIMITS.items())),
         "指数门上限（log(秒) vs log(规模) 最小二乘斜率）"),
    ]
    lines = ["| 常量（`tensorcad.perf_benchmark`） | 值 | 作用 |", "|---|---|---|"]
    for name, value, why in facts:
        lines.append(f"| {name} | {value} | {why} |")
    return "\n".join(lines)


#: table name -> (artifact keys it consumes, renderer)
TABLE_SPECS: dict[str, tuple[tuple[str, ...], Any]] = {
    "mesh-table": (("mesh",), render_mesh_table),
    "voxel-table": (("voxel",), render_voxel_table),
    "budget-table": (("mesh", "voxel"), render_budget_table),
    "fleet3pt-table": (("fleet3pt",), render_fleet3pt_table),
    "fleet3pt-gates": (("fleet3pt",), render_fleet3pt_gates),
}

#: Table names the document can embed, in display order.
TABLE_NAMES: tuple[str, ...] = tuple(TABLE_SPECS)

#: Tables of ``docs/perf_benchmark.md`` (the gate-conventions manual, not the
#: measurement report).  Kept in a separate registry because a document must not
#: be required to embed tables it does not need: an absent marker counts as
#: drift, and mixing the two bundles would make both documents fail.
BENCHMARK_DOC_SPECS: dict[str, tuple[tuple[str, ...], Any]] = {
    "baseline-gate-table": (("baseline",), render_baseline_gate_table),
    "code-facts-table": ((), render_code_facts),
}

#: Table names ``docs/perf_benchmark.md`` must embed, in display order.
BENCHMARK_DOC_TABLE_NAMES: tuple[str, ...] = tuple(BENCHMARK_DOC_SPECS)

#: Every renderer, from both registries.
ALL_TABLE_SPECS: dict[str, tuple[tuple[str, ...], Any]] = {
    **TABLE_SPECS,
    **BENCHMARK_DOC_SPECS,
}


def render_block(name: str, artifacts: Mapping[str, Mapping[str, Any]]) -> str:
    """Rendered body (markers excluded) of table ``name``."""
    try:
        keys, renderer = ALL_TABLE_SPECS[name]
    except KeyError:  # pragma: no cover - guard for typos
        raise ValueError(f"unknown perf-evidence table {name!r}") from None
    missing = [k for k in keys if k not in artifacts]
    if missing:
        raise ValueError(f"{name} needs artifact(s) {missing}")
    return renderer(*[artifacts[k] for k in keys])


def extract_block(doc_text: str, name: str) -> str | None:
    """Body between the markers for ``name``, or ``None`` when absent."""
    begin, end = marker_begin(name), marker_end(name)
    i = doc_text.find(begin)
    if i < 0:
        return None
    start = i + len(begin)
    j = doc_text.find(end, start)
    if j < 0:
        return None
    return doc_text[start:j].strip("\n")


def replace_block(doc_text: str, name: str, body: str) -> str:
    """Return ``doc_text`` with the ``name`` block body replaced by ``body``."""
    begin, end = marker_begin(name), marker_end(name)
    i = doc_text.find(begin)
    if i < 0:
        raise ValueError(f"no begin marker for {name!r} in document")
    start = i + len(begin)
    j = doc_text.find(end, start)
    if j < 0:
        raise ValueError(f"no end marker for {name!r} in document")
    return doc_text[:start] + "\n" + body.rstrip("\n") + "\n" + doc_text[j:]


def check_document(
    doc_text: str,
    artifacts: Mapping[str, Mapping[str, Any]],
    names: Sequence[str] | None = None,
) -> list[str]:
    """Return the drift problems between ``doc_text`` and ``artifacts``.

    An empty list means every marked table is exactly what the committed
    measurement JSON renders to.  A missing marker is reported too: a table that
    silently stopped being rendered is the failure mode this module prevents.
    """
    problems: list[str] = []
    for name in TABLE_NAMES if names is None else names:
        body = extract_block(doc_text, name)
        if body is None:
            problems.append(f"{name}: markers not found in document")
            continue
        try:
            expected = render_block(name, artifacts)
        except ValueError as exc:
            problems.append(f"{name}: {exc}")
            continue
        if body.strip() != expected.strip():
            problems.append(
                f"{name}: document body differs from rendered artifact "
                f"({len(body.splitlines())} vs {len(expected.splitlines())} lines)"
            )
    return problems


def load_artifact(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def render_document(
    doc_text: str,
    artifacts: Mapping[str, Mapping[str, Any]],
    names: Sequence[str] | None = None,
) -> str:
    """Return ``doc_text`` with every marked block re-rendered from artifacts."""
    out = doc_text
    for name in TABLE_NAMES if names is None else names:
        keys = ALL_TABLE_SPECS[name][0]
        if all(k in artifacts for k in keys):
            out = replace_block(out, name, render_block(name, artifacts))
    return out