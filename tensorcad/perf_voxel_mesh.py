"""Voxelisation *from an explicit triangle mesh*: the 10^6--10^7 triangle regime.

Why this is a separate ladder
-----------------------------
The ``voxelize`` sub-task in :mod:`tensorcad.perf_benchmark` samples
``SuboffHullField`` -- an analytic field. Its cost is closed-form per point, so
its ladder is a *grid-resolution* one and its triangle column is structurally
zero (0 triangles involved). That answers "how long does it take to sample an
analytic SDF onto 10^8 cells", which is **not** the question a solver pipeline
asks: an imported CAD mesh, a repaired STL or a B-rep tessellation at 10^6--10^7
triangles has no closed form, and the SDF has to be derived from the triangles.

This module measures that second regime. Every rung builds a
:class:`~tensorcad.mesh_field.MeshSDFField` over an explicit mesh (one connected
closed body, from ``generate_hull`` or from 1-to-4 midpoint subdivision of a
coarser hull) and voxelises it with an exact Euclidean SDF truncated to a band of
``band_factor`` voxel widths -- the narrow-band rule real voxelisers use, which
is exact wherever it matters for the solid fraction and caps the cost of the
deep-interior cells that would otherwise have to clear a half-metre ball of
triangles each.

Two size axes are swept, because the cost has two factors and only one of them
was ever measured before:

``triangles``  mesh size at a fixed grid (the ladder that reaches 10^7);
``cells``      grid size at a fixed mesh.

Gates
-----
* ``scaling_voxelize_mesh`` -- complexity exponent against triangle count at a
  fixed grid (limit 1.5, i.e. anything close to quadratic is RED).
* ``scaling_voxelize_mesh_cells`` -- exponent against cell count at a fixed mesh
  (limit 1.5; the analytic path measured 0.98, so ~1 is the honest expectation).
* ``voxel_volume_accuracy`` -- voxelised volume vs the analytic hull volume,
  flagged above 3 % relative error. A fast but wrong voxelisation cannot pass.

Nothing here is extrapolated: rungs above ``max_triangles`` are recorded
``not_run`` with the expected triangle count, and a rung that exceeds the
per-rung budget is reported as a measured row plus a *cluster* warning with the
number that produced it.

Usage
-----
::

    # mechanism smoke (seconds)
    python scripts/perf_voxelize_mesh.py --bases 41:13 --levels 0,1 --grids 8:4:2,16:8:4 \\
        --out /tmp/voxel_mesh.json

    # the real run (minutes; submit through the compute scheduler)
    python scripts/perf_voxelize_mesh.py --bases 41:13,401:129 --levels 0,1,2,3 \\
        --grids 16:8:4,24:12:6,32:16:8,48:24:12 \\
        --output benchmarks/perf_scale/artifacts/ladder_voxelize_mesh.json
"""
from __future__ import annotations

import math
import platform
import sys
import time
from typing import Any, Iterable, Sequence

from .geometry import Mesh
from .hull import HullParameters, HullModel, analytic_volume, generate_hull
from .mesh_field import MeshSDFField
from .mesh_subdivision import expected_face_count, subdivide_mesh
from .perf_benchmark import (
    live_rss_mb,
    peak_rss_mb,
    _environment,
)
from .spatial_index import TriangleGrid

SCHEMA = "tensorcad.perf_voxel_mesh/1"

#: One coarse base (41x13 = 6.6e3 triangles) is enough to reach 10^7 by
#: subdivision: 6.6e3 -> 2.6e4 -> 1.0e5 -> 4.1e5 -> 1.6e6 -> 6.6e6 -> 2.6e7.
DEFAULT_BASES: tuple[tuple[int, int], ...] = ((41, 13), (401, 129))
DEFAULT_LEVELS: tuple[int, ...] = (0, 1, 2, 3)
#: Grids as (nx, ny, nz). Anisotropic, because a 10 m x 2 m x 0.8 m hull in a
#: cubic grid wastes 4/5 of the cells and still under-resolves the length.
DEFAULT_GRIDS: tuple[tuple[int, int, int], ...] = (
    (8, 4, 2), (16, 8, 4), (24, 12, 6), (32, 16, 8), (48, 24, 12),
)
DEFAULT_HULL = HullParameters(length=10.0, beam=2.0, draft=0.8)

#: Exponent limits, matching the mesh ladders (a quadratic path is a regression).
SCALING_LIMIT_TRIANGLES = 1.5
SCALING_LIMIT_CELLS = 1.5
#: Fewest distinct sizes needed before an exponent is fitted at all.
MIN_SCALING_POINTS = 3
#: Per-rung wall clock above which the row says "cluster", not "one machine".
DEFAULT_BUDGET_SECONDS = 600.0
#: Volume error above which a rung is RED.
VOLUME_TOLERANCE = 0.03


def _to_mesh(obj: HullModel | Mesh) -> Mesh:
    return obj.mesh if isinstance(obj, HullModel) else obj


def _grid_spacing(mesh: Mesh, grid: Sequence[int]) -> tuple[float, float, float]:
    vertices = mesh.vertices
    lo = [min(v[axis] for v in vertices) for axis in range(3)]
    hi = [max(v[axis] for v in vertices) for axis in range(3)]
    return tuple(
        max(hi[axis] - lo[axis], 1e-12) / int(grid[axis]) for axis in range(3)
    )  # type: ignore[return-value]


def _row(**fields: Any) -> dict[str, Any]:
    row = {
        "subtask": "voxelize_mesh",
        "nx": 0,
        "nz": 0,
        "levels": 0,
        "grid": [0, 0, 0],
        "grid_cells": 0,
        "triangles": 0,
        "vertices": 0,
        "seconds": None,
        "status": "measured",
        "queries_per_second": None,
        "us_per_query": None,
        "band": None,
        "index_cells": 0,
        "index_entries": 0,
        "index_bytes": 0,
        "build_seconds": None,
        "peak_rss_mb": None,
        "rss_mb": None,
        "mesh_mb": None,
        "voxel_volume": None,
        "analytic_volume": None,
        "volume_rel_err": None,
        "note": "",
    }
    row.update(fields)
    return row


def _voxelize_rung(
    mesh: Mesh,
    grid: tuple[int, int, int],
    analytic: float,
    *,
    band_factor: float,
    nx: int,
    nz: int,
    levels: int,
    closed_checks: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Voxelise one mesh at one grid, timing query path and index build apart."""
    from .voxel import to_voxel_grid

    spacing = _grid_spacing(mesh, grid)
    band = band_factor * max(spacing)

    t0 = time.perf_counter()
    index = TriangleGrid(mesh)
    build_seconds = time.perf_counter() - t0
    stats = index.stats()
    field = MeshSDFField(mesh, grid=index, band=band)

    cells = int(grid[0]) * int(grid[1]) * int(grid[2])
    t0 = time.perf_counter()
    voxels = to_voxel_grid(field, grid, bounds=None)
    seconds = time.perf_counter() - t0

    volume = voxels.volume()
    rel_err = (volume - analytic) / analytic if analytic else float("nan")

    row = _row(
        nx=nx,
        nz=nz,
        levels=levels,
        grid=[int(v) for v in grid],
        grid_cells=cells,
        triangles=len(mesh.faces),
        vertices=len(mesh.vertices),
        seconds=round(seconds, 6),
        queries_per_second=round(cells / seconds, 1) if seconds > 0 else None,
        us_per_query=round(1e6 * seconds / cells, 1) if cells else None,
        band=round(band, 6),
        index_cells=int(stats["cells"]),
        index_entries=int(stats["entries"]),
        index_bytes=int(stats["index_bytes"]),
        build_seconds=round(build_seconds, 6),
        peak_rss_mb=peak_rss_mb(),
        rss_mb=live_rss_mb(),
        mesh_mb=sum(len(v) * 8 for v in mesh.vertices) / 1e6
        + len(mesh.faces) * 12 / 1e6,
        voxel_volume=round(volume, 9),
        analytic_volume=round(analytic, 9),
        volume_rel_err=round(rel_err, 9) if analytic else None,
        note=(
            f"exact-Euclidean SDF truncated at band={band:.4f} m "
            f"({band_factor:g} voxel widths); band >= half a cell so the solid "
            f"fraction is unaffected"
        ),
    )
    if closed_checks:
        row.update(closed_checks)
    return row


def _loglog_slope(points: list[tuple[float, float]]) -> tuple[float, float]:
    xs = [math.log(x) for x, _ in points]
    ys = [math.log(y) for _, y in points]
    n = len(points)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0.0:
        raise ValueError("degenerate size axis")
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    intercept = my - slope * mx
    ss_res = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - my) ** 2 for y in ys)
    return slope, (1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0)


def scaling_gates(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """Fit the two complexity exponents and gate the volume accuracy per rung."""
    gates: list[dict[str, Any]] = []

    # Triangle axis: group by grid, keep grids with >= 3 measured rungs.
    by_grid: dict[tuple[int, ...], list[dict[str, Any]]] = {}
    for row in rows:
        if row["status"] != "measured" or not row["seconds"]:
            continue
        by_grid.setdefault(tuple(row["grid"]), []).append(row)
    tri_fits: list[tuple[float, float, int, str]] = []
    for grid, group in sorted(by_grid.items()):
        sizes = sorted({r["triangles"] for r in group})
        if len(sizes) < MIN_SCALING_POINTS:
            continue
        points = [(float(r["triangles"]), float(r["seconds"])) for r in group]
        sizes_label = f"{sizes[0]:,}..{sizes[-1]:,}"
        slope, r2 = _loglog_slope(points)
        tri_fits.append((slope, r2, len(sizes), f"grid={grid[0]}x{grid[1]}x{grid[2]} sizes {sizes_label}"))
    if tri_fits:
        worst = max(tri_fits, key=lambda f: f[0])
        gates.append({
            "name": "scaling_voxelize_mesh",
            "status": "measured",
            "passed": worst[0] <= SCALING_LIMIT_TRIANGLES,
            "detail": (
                f"worst exponent={worst[0]:.2f} (limit {SCALING_LIMIT_TRIANGLES:.2f}), "
                f"R^2={worst[1]:.3f}, {worst[2]} sizes, {worst[3]}"
            ),
            "ratio": round(worst[0], 3),
        })
    else:
        gates.append({
            "name": "scaling_voxelize_mesh",
            "status": "not_run",
            "passed": None,
            "detail": (
                f"not_run: need >= {MIN_SCALING_POINTS} measured triangle sizes at one "
                f"grid to fit an exponent"
            ),
            "ratio": None,
        })

    # Cell axis: group by mesh (triangles), keep groups with >= 3 measured rungs.
    by_mesh: dict[int, list[dict[str, Any]]] = {}
    for row in rows:
        if row["status"] != "measured" or not row["seconds"]:
            continue
        by_mesh.setdefault(int(row["triangles"]), []).append(row)
    cell_fits: list[tuple[float, float, int, str]] = []
    for triangles, group in sorted(by_mesh.items()):
        sizes = sorted({r["grid_cells"] for r in group})
        if len(sizes) < MIN_SCALING_POINTS:
            continue
        points = [(float(r["grid_cells"]), float(r["seconds"])) for r in group]
        slope, r2 = _loglog_slope(points)
        cell_fits.append((slope, r2, len(sizes), f"mesh {triangles:,} tri, cells {sizes[0]:,}..{sizes[-1]:,}"))
    if cell_fits:
        worst = max(cell_fits, key=lambda f: f[0])
        gates.append({
            "name": "scaling_voxelize_mesh_cells",
            "status": "measured",
            "passed": worst[0] <= SCALING_LIMIT_CELLS,
            "detail": (
                f"worst exponent={worst[0]:.2f} (limit {SCALING_LIMIT_CELLS:.2f}), "
                f"R^2={worst[1]:.3f}, {worst[2]} sizes, {worst[3]}"
            ),
            "ratio": round(worst[0], 3),
        })
    else:
        gates.append({
            "name": "scaling_voxelize_mesh_cells",
            "status": "not_run",
            "passed": None,
            "detail": (
                f"not_run: need >= {MIN_SCALING_POINTS} measured cell counts at one mesh"
            ),
            "ratio": None,
        })

    # Accuracy: every measured rung with an analytic reference must be < 3 %.
    checked = [
        r for r in rows
        if r["status"] == "measured" and r["volume_rel_err"] is not None
    ]
    if checked:
        worst = max(checked, key=lambda r: abs(r["volume_rel_err"]))
        gates.append({
            "name": "voxel_volume_accuracy",
            "status": "measured",
            "passed": abs(worst["volume_rel_err"]) <= VOLUME_TOLERANCE,
            "detail": (
                f"worst |rel err|={abs(worst['volume_rel_err']) * 100:.3f}% "
                f"(tolerance {VOLUME_TOLERANCE * 100:.0f}%) at "
                f"{worst['triangles']:,} tri, grid {tuple(worst['grid'])}"
            ),
            "ratio": round(abs(worst["volume_rel_err"]), 6),
        })
    else:
        gates.append({
            "name": "voxel_volume_accuracy",
            "status": "not_run",
            "passed": None,
            "detail": "not_run: no measured rung carries an analytic reference volume",
            "ratio": None,
        })
    return gates


def run_voxelize_mesh_benchmark(
    bases: Iterable[Sequence[int]] = DEFAULT_BASES,
    levels_ladder: Iterable[int] = DEFAULT_LEVELS,
    grids: Iterable[Sequence[int]] = DEFAULT_GRIDS,
    hull: HullParameters = DEFAULT_HULL,
    band_factor: float = 1.0,
    max_triangles: int | None = 30_000_000,
    budget_seconds: float = DEFAULT_BUDGET_SECONDS,
    command: str | None = None,
) -> dict[str, Any]:
    """Sweep both size axes and return a JSON-ready evidence record."""
    if band_factor <= 0:
        raise ValueError("band_factor must be positive")
    levels_sorted = sorted({int(level) for level in levels_ladder})
    if levels_sorted and levels_sorted[0] < 0:
        raise ValueError("levels must be >= 0")
    grid_list = [tuple(int(v) for v in g) for g in grids]
    for g in grid_list:
        if len(g) != 3 or min(g) < 1:
            raise ValueError(f"grids must be positive (nx, ny, nz) triples, got {g!r}")

    rows: list[dict[str, Any]] = []
    cluster_warnings: list[str] = []
    t_start = time.perf_counter()
    analytic = analytic_volume(
        hull.length, hull.beam, hull.draft,
        fore_exponent=hull.fore_exponent, aft_exponent=hull.aft_exponent,
    )

    for nx, nz in bases:
        nx, nz = int(nx), int(nz)
        base = _to_mesh(generate_hull(HullParameters(
            length=hull.length, beam=hull.beam, draft=hull.draft, nx=nx, nz=nz,
        )))
        base_triangles = len(base.faces)
        mesh = base
        for level in levels_sorted:
            if level:
                mesh = subdivide_mesh(mesh, 1)
            triangles = len(mesh.faces)
            if max_triangles is not None and triangles > max_triangles:
                rows.append(_row(
                    nx=nx, nz=nz, levels=level, triangles=triangles,
                    status="not_run",
                    note=(
                        f"{triangles:,} triangles exceeds max_triangles="
                        f"{max_triangles:,}; not voxelised"
                    ),
                ))
                continue
            for grid in grid_list:
                try:
                    row = _voxelize_rung(
                        mesh, grid, analytic,
                        band_factor=band_factor, nx=nx, nz=nz, levels=level,
                    )
                except Exception as exc:  # noqa: BLE001 - recorded, never swallowed
                    rows.append(_row(
                        nx=nx, nz=nz, levels=level, grid=list(grid), triangles=triangles,
                        status="failed", note=f"{type(exc).__name__}: {exc}",
                    ))
                    continue
                rows.append(row)
                if row["seconds"] and row["seconds"] > budget_seconds:
                    cluster_warnings.append(
                        f"voxelize_mesh at {triangles:,} tri / grid "
                        f"{grid[0]}x{grid[1]}x{grid[2]} ({row['grid_cells']:,} cells) took "
                        f"{row['seconds']:.1f}s > {budget_seconds:.0f}s budget — "
                        f"submit larger rungs to the cluster"
                    )
                if abs(row["volume_rel_err"] or 0.0) > VOLUME_TOLERANCE:
                    cluster_warnings.append(
                        f"voxelize_mesh volume error {row['volume_rel_err'] * 100:.2f}% "
                        f"at {triangles:,} tri / grid {grid[0]}x{grid[1]}x{grid[2]} "
                        f"exceeds {VOLUME_TOLERANCE * 100:.0f}%"
                    )

    gates = scaling_gates(rows)
    elapsed = time.perf_counter() - t_start
    return {
        "schema": SCHEMA,
        "passed": all(g["passed"] for g in gates if g["status"] == "measured"),
        "verdict": "pass" if all(
            g["passed"] for g in gates if g["status"] == "measured"
        ) else "fail",
        "gates_checked": sum(1 for g in gates if g["status"] == "measured"),
        "gates_not_run": sum(1 for g in gates if g["status"] == "not_run"),
        "seconds": round(elapsed, 3),
        "command": command or "",
        "band_factor": band_factor,
        "grids": [list(g) for g in grid_list],
        "levels": levels_sorted,
        "bases": [[int(nx), int(nz)] for nx, nz in bases],
        "solver": {
            "sdf": "exact Euclidean, truncated to band",
            "sign": "ray-crossing parity",
            "field": "tensorcad.mesh_field.MeshSDFField",
        },
        "environment": _environment(),
        "cluster_warnings": cluster_warnings,
        "gates": gates,
        "results": rows,
    }


# ---------------------------------------------------------------------------
# Markdown rendering (the doc's table is generated from the artifact, not typed)
# ---------------------------------------------------------------------------

MESH_TABLE_BEGIN = "<!-- perf-meshvoxel:tri-table begin -->"
MESH_TABLE_END = "<!-- perf-meshvoxel:tri-table end -->"
GRID_TABLE_BEGIN = "<!-- perf-meshvoxel:grid-table begin -->"
GRID_TABLE_END = "<!-- perf-meshvoxel:grid-table end -->"
GATE_TABLE_BEGIN = "<!-- perf-meshvoxel:gates begin -->"
GATE_TABLE_END = "<!-- perf-meshvoxel:gates end -->"


def _table(header: Sequence[str], align: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join(align) + "|"]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def render_triangle_table(artifact: dict[str, Any], grid: Sequence[int] | None = None) -> str:
    """Rows of the triangle sweep at one fixed grid (default: the finest with rows)."""
    rows = [r for r in artifact["results"] if r["status"] == "measured"]
    if not rows:
        return "_(no measured rungs)_"
    if grid is None:
        # The grid with the most measured rungs: a table with a single row cannot
        # show a trend, so "the finest grid" is the wrong default.
        counts: dict[tuple[int, ...], int] = {}
        for row in rows:
            key = tuple(row["grid"])
            counts[key] = counts.get(key, 0) + 1
        grid = max(counts, key=lambda key: (counts[key], key))
    grid = tuple(grid)
    body = []
    for row in sorted(
        (r for r in rows if tuple(r["grid"]) == grid), key=lambda r: r["triangles"]
    ):
        body.append([
            f"{row['triangles']:,}",
            f"{row['seconds']:.3f}",
            f"{row['us_per_query']:.0f}",
            f"{row['band']:.3f}",
            f"{row['index_cells']:,}",
            f"{row['index_bytes'] / 2**20:.1f}",
            f"{row['volume_rel_err'] * 100:+.3f}%",
        ])
    return _table(
        ["三角面", "秒", "µs/查询", "band m", "索引 cells", "索引 MiB", "体素体积误差"],
        ["---:", "---:", "---:", "---:", "---:", "---:", "---:"],
        body,
    )


def render_cell_table(artifact: dict[str, Any], triangles: int | None = None) -> str:
    """Rows of the cell sweep at one fixed mesh (default: the largest with rows)."""
    rows = [r for r in artifact["results"] if r["status"] == "measured"]
    if not rows:
        return "_(no measured rungs)_"
    if triangles is None:
        counts: dict[int, int] = {}
        for row in rows:
            counts[row["triangles"]] = counts.get(row["triangles"], 0) + 1
        triangles = max(counts, key=lambda key: (counts[key], key))
    body = []
    for row in sorted(
        (r for r in rows if r["triangles"] == triangles), key=lambda r: r["grid_cells"]
    ):
        body.append([
            f"{row['grid'][0]}×{row['grid'][1]}×{row['grid'][2]}",
            f"{row['grid_cells']:,}",
            f"{row['seconds']:.3f}",
            f"{row['us_per_query']:.0f}",
            f"{row['band']:.3f}",
            f"{row['volume_rel_err'] * 100:+.3f}%",
        ])
    return _table(
        ["网格", "体素", "秒", "µs/查询", "band m", "体素体积误差"],
        ["---:", "---:", "---:", "---:", "---:", "---:"],
        body,
    )


def render_gate_table(artifact: dict[str, Any]) -> str:
    body = []
    for gate in artifact["gates"]:
        verdict = (
            "pass" if gate["passed"] else "RED"
        ) if gate["status"] == "measured" else "not_run"
        body.append([f"`{gate['name']}`", gate["status"], verdict, gate["detail"]])
    return _table(["门", "状态", "判定", "证据"], ["---", "---", "---", "---"], body)


def render_section(artifact: dict[str, Any]) -> str:
    """The full (marker-delimited) Markdown block the doc embeds."""
    parts = [
        f"`{artifact['schema']}` · {artifact['verdict']} · "
        f"gates checked={artifact['gates_checked']} not_run={artifact['gates_not_run']} · "
        f"wall {artifact['seconds'] / 60:.1f} min · band={artifact['band_factor']:g} 体宽",
        "",
        "#### 三角面轴（固定网格，窄带 1 体宽）",
        "",
        MESH_TABLE_BEGIN,
        render_triangle_table(artifact),
        MESH_TABLE_END,
        "",
        "#### 体素轴（固定网格数 vs 固定网格分辨率，最大单体网格）",
        "",
        GRID_TABLE_BEGIN,
        render_cell_table(artifact),
        GRID_TABLE_END,
        "",
        "#### 门（`not_run` 不是 pass）",
        "",
        GATE_TABLE_BEGIN,
        render_gate_table(artifact),
        GATE_TABLE_END,
    ]
    return "\n".join(parts)


def check_section(doc_text: str, artifact: dict[str, Any]) -> list[str]:
    """Compare the doc's rendered blocks against the artifact; return mismatches."""
    problems: list[str] = []
    blocks = [
        (MESH_TABLE_BEGIN, MESH_TABLE_END, render_triangle_table(artifact)),
        (GRID_TABLE_BEGIN, GRID_TABLE_END, render_cell_table(artifact)),
        (GATE_TABLE_BEGIN, GATE_TABLE_END, render_gate_table(artifact)),
    ]
    for begin, end, expected in blocks:
        if begin not in doc_text or end not in doc_text:
            problems.append(f"missing marker block {begin}")
            continue
        actual = doc_text.split(begin, 1)[1].split(end, 1)[0].strip()
        if actual != expected.strip():
            problems.append(f"table drift in {begin}")
    return problems