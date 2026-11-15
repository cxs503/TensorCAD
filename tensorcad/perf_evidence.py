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
from typing import Any, Mapping

#: Measured sub-tasks, in display order, with the axis they scale in.
MESH_ORDER = ("load", "point_locate", "point_locate_accel", "volume")
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


#: table name -> (artifact keys it consumes, renderer)
TABLE_SPECS: dict[str, tuple[tuple[str, ...], Any]] = {
    "mesh-table": (("mesh",), render_mesh_table),
    "voxel-table": (("voxel",), render_voxel_table),
    "budget-table": (("mesh", "voxel"), render_budget_table),
}

#: Table names the document can embed, in display order.
TABLE_NAMES: tuple[str, ...] = tuple(TABLE_SPECS)


def render_block(name: str, artifacts: Mapping[str, Mapping[str, Any]]) -> str:
    """Rendered body (markers excluded) of table ``name``."""
    try:
        keys, renderer = TABLE_SPECS[name]
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
) -> list[str]:
    """Return the drift problems between ``doc_text`` and ``artifacts``.

    An empty list means every marked table is exactly what the committed
    measurement JSON renders to.  A missing marker is reported too: a table that
    silently stopped being rendered is the failure mode this module prevents.
    """
    problems: list[str] = []
    for name in TABLE_NAMES:
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
) -> str:
    """Return ``doc_text`` with every marked block re-rendered from artifacts."""
    out = doc_text
    for name in TABLE_NAMES:
        keys = TABLE_SPECS[name][0]
        if all(k in artifacts for k in keys):
            out = replace_block(out, name, render_block(name, artifacts))
    return out