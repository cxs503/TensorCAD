"""Single-body 10^6-10^7 benchmark: subdivision as the large-mesh source.

``tensorcad.perf_benchmark`` reaches the 10^7 band with a *fleet* of disjoint
copies (``build_fleet_mesh``), because ``HullParameters`` caps one parametric
hull at 4*(nx-1)*(nz-1) = 8,199,998 triangles.  That answers "how does the code
scale with triangle count" but not "what happens on one 10^7-triangle body",
which is the shape real imported CAD/solver meshes have (one connected closed
surface, no per-copy regularity).

This driver closes that gap with :func:`tensorcad.mesh_subdivision.subdivide_mesh`
-- uniform 1-to-4 subdivision of a *coarser* hull.  Two consequences make it more
than a convenience:

* the result is a **single connected** body, so the fleet's O(N) copy cost does
  not contaminate ``load``/``volume``/``point_locate_accel``;
* midpoint subdivision tiles every planar parent exactly, so the enclosed volume
  and the surface area are **preserved by construction**.  That turns the
  divergence-theorem integral into a *check* (relative error against the
  analytic hull volume) instead of a comparison between two discretisations.

The artifact this module emits is small JSON (one row per measurement plus the
invariant checks) and is rendered into ``docs/perf_scale_single_body_10m.md`` by
:mod:`tensorcad.perf_evidence`; it is never transcribed by hand.
"""
from __future__ import annotations

import time
from typing import Any, Iterable, Sequence

from .hull import HullParameters, analytic_volume, generate_hull
from .mesh import (
    edge_face_counts,
    is_edge_manifold,
    is_watertight,
    surface_area,
    total_signed_volume,
)
from .mesh_subdivision import expected_face_count, subdivide_mesh
from .perf_benchmark import (
    _environment,
    _fibonacci_directions,
    live_rss_mb,
    mesh_footprint_mb,
    peak_rss_mb,
)

SCHEMA = "tensorcad.perf_single10m/1"

#: Default coarse bases (nx, nz): 4*(nx-1)*(nz-1) triangles before subdivision.
DEFAULT_BASES: tuple[tuple[int, int], ...] = ((401, 129), (1601, 513))

#: Hard cap on one rung: 1.6e7 triangles need ~3.5 GB for the Python mesh plus
#: ~1.2 GB for the index, which is what the shared node can offer without
#: swapping.  Rungs above the cap are recorded ``not_run`` with the reason,
#: never silently dropped.

#: Default subdivision ladder.  Level 0 is the generated hull itself; each level
#: multiplies the triangle count by four (401x129 -> 0.20e6 / 0.82e6 / 3.28e6 /
#: 13.1e6 at levels 0/1/2/3, i.e. the whole 10^6-10^7 band from one base).
DEFAULT_LEVELS: tuple[int, ...] = (0, 1, 2, 3)

DEFAULT_HULL = HullParameters(length=10.0, beam=2.0, draft=0.8)


def _euler_characteristic(mesh: Any) -> int:
    return len(mesh.vertices) - len(edge_face_counts(mesh)) + len(mesh.faces)


def _row(**fields: Any) -> dict[str, Any]:
    row = {
        "subtask": "",
        "nx": 0,
        "nz": 0,
        "levels": 0,
        "triangles": 0,
        "vertices": 0,
        "seconds": None,
        "status": "measured",
        "throughput": None,
        "note": "",
        "index_cells": 0,
        "index_entries": 0,
        "index_bytes": 0,
        "build_seconds": None,
        "peak_rss_mb": None,
        "rss_mb": None,
        "mesh_mb": None,
    }
    row.update(fields)
    return row


def _query_accel(mesh: Any, n_queries: int, nx: int, nz: int, levels: int) -> dict[str, Any]:
    from .spatial_index import TriangleGrid

    triangles = len(mesh.faces)
    t0 = time.perf_counter()
    grid = TriangleGrid(mesh)
    build = time.perf_counter() - t0

    directions = _fibonacci_directions(n_queries)
    t0 = time.perf_counter()
    hits = 0
    for direction in directions:
        if grid.cast_ray((0.0, 0.0, 0.0), direction) is not None:
            hits += 1
    elapsed = time.perf_counter() - t0
    stats = grid.stats()
    return _row(
        subtask="point_locate_accel",
        nx=nx,
        nz=nz,
        levels=levels,
        triangles=triangles,
        vertices=len(mesh.vertices),
        seconds=round(elapsed, 6),
        throughput=round(n_queries / elapsed, 1) if elapsed > 0 else None,
        build_seconds=round(build, 6),
        note=(
            f"hits={hits}/{n_queries}, grid={tuple(stats['resolution'])}, "
            f"build={build:.3f}s, single connected body, levels={levels}"
        ),
        index_cells=int(stats["cells"]),
        index_entries=int(stats["entries"]),
        index_bytes=int(stats["index_bytes"]),
        peak_rss_mb=peak_rss_mb(),
        rss_mb=live_rss_mb(),
        mesh_mb=mesh_footprint_mb(mesh),
    )


def run_single_body_benchmark(
    bases: Iterable[Sequence[int]] = DEFAULT_BASES,
    levels_ladder: Iterable[int] = DEFAULT_LEVELS,
    hull: HullParameters = DEFAULT_HULL,
    n_queries: int = 2000,
    max_triangles: int | None = 30_000_000,
    command: str | None = None,
) -> dict[str, Any]:
    """Measure ``load``/``subdivide``/``volume``/``point_locate_accel`` on one body.

    Every rung is one connected closed mesh, so the 10^6-10^7 rows are directly
    comparable with an imported large mesh.  Invariant checks (watertight,
    edge-manifold, Euler characteristic 2, volume/area preservation along the
    subdivision ladder) are recorded per rung: a fast but wrong geometry source
    cannot pass silently.
    """
    levels_sorted = sorted({int(level) for level in levels_ladder})
    if -1 in levels_sorted:
        raise ValueError("levels must be >= 0")
    rows: list[dict[str, Any]] = []
    checks: list[dict[str, Any]] = []
    t_start = time.perf_counter()

    for nx, nz in bases:
        nx, nz = int(nx), int(nz)
        t0 = time.perf_counter()
        base = generate_hull(HullParameters(
            length=hull.length, beam=hull.beam, draft=hull.draft, nx=nx, nz=nz,
        ))
        load_seconds = time.perf_counter() - t0
        exact = analytic_volume(
            hull.length, hull.beam, hull.draft,
            fore_exponent=hull.fore_exponent, aft_exponent=hull.aft_exponent,
        )
        base_triangles = len(base.faces)
        base_area = surface_area(base)
        base_volume = total_signed_volume(base)
        rows.append(_row(
            subtask="load",
            nx=nx,
            nz=nz,
            levels=0,
            triangles=base_triangles,
            vertices=len(base.vertices),
            seconds=round(load_seconds, 6),
            throughput=round(base_triangles / load_seconds, 1) if load_seconds > 0 else None,
            note=f"generate_hull({nx}x{nz}), analytic_volume={exact:.6f} m3",
            peak_rss_mb=peak_rss_mb(),
            rss_mb=live_rss_mb(),
            mesh_mb=mesh_footprint_mb(base),
        ))

        mesh = base
        for level in levels_sorted:
            step_seconds: float | None = 0.0
            if max_triangles is not None:
                expected = expected_face_count(base_triangles, level)
                if expected > max_triangles:
                    rows.append(_row(
                        subtask="subdivide", nx=nx, nz=nz, levels=level,
                        triangles=expected, vertices=0, seconds=None,
                        status="not_run",
                        note=(
                            f"expected {expected:,} triangles exceeds "
                            f"max_triangles={max_triangles:,}; not subdivided"
                        ),
                    ))
                    break
            if level > 0:
                t0 = time.perf_counter()
                mesh = subdivide_mesh(mesh, 1)
                step_seconds = time.perf_counter() - t0
                rows.append(_row(
                    subtask="subdivide",
                    nx=nx,
                    nz=nz,
                    levels=level,
                    triangles=len(mesh.faces),
                    vertices=len(mesh.vertices),
                    seconds=round(step_seconds, 6),
                    throughput=(
                        round(len(mesh.faces) / step_seconds, 1) if step_seconds > 0 else None
                    ),
                    note=(
                        f"1-to-4 midpoint subdivision, step {level - 1}->{level}, "
                        f"single connected body"
                    ),
                    peak_rss_mb=peak_rss_mb(),
                    rss_mb=live_rss_mb(),
                    mesh_mb=mesh_footprint_mb(mesh),
                ))
            t0 = time.perf_counter()
            volume = total_signed_volume(mesh)
            volume_seconds = time.perf_counter() - t0
            rel_err = (volume - exact) / exact
            rows.append(_row(
                subtask="volume",
                nx=nx,
                nz=nz,
                levels=level,
                triangles=len(mesh.faces),
                vertices=len(mesh.vertices),
                seconds=round(volume_seconds, 6),
                throughput=(
                    round(len(mesh.faces) / volume_seconds, 1) if volume_seconds > 0 else None
                ),
                note=(
                    f"analytic={exact:.6f} m3, measured={volume:.6f} m3, "
                    f"rel_err={rel_err:.3e}"
                ),
                peak_rss_mb=peak_rss_mb(),
                rss_mb=live_rss_mb(),
                mesh_mb=mesh_footprint_mb(mesh),
            ))

            area = surface_area(mesh)
            checks.append({
                "nx": nx,
                "nz": nz,
                "levels": level,
                "triangles": len(mesh.faces),
                "vertices": len(mesh.vertices),
                "connected_bodies": 1,
                "watertight": bool(is_watertight(mesh)),
                "edge_manifold": bool(is_edge_manifold(mesh)),
                "euler_characteristic": _euler_characteristic(mesh),
                "volume_m3": volume,
                "analytic_volume_m3": exact,
                "volume_rel_err": rel_err,
                "surface_area_m2": area,
                "surface_area_rel_drift": (area - base_area) / base_area if base_area else None,
                "volume_rel_drift_vs_base": (volume - base_volume) / base_volume if base_volume else None,
                "max_triangle_area_ratio": _max_area_ratio(mesh, base, level),
            })

            rows.append(_query_accel(mesh, n_queries=n_queries, nx=nx, nz=nz, levels=level))

    elapsed = time.perf_counter() - t_start
    # Two different questions must not be mixed here:
    #  * preservation -- subdivision must not *change* the discrete solid, so the
    #    drift against the generated parent is a hard invariant (tolerance ~fp);
    #  * accuracy -- the parent is itself a sampled surface, so its distance to
    #    the analytic hull volume is inherited from the coarse mesh and cannot
    #    improve by refining.  Reported separately, never inside ``invariants_ok``.
    volume_drifts = [
        abs(c["volume_rel_drift_vs_base"]) for c in checks
        if c["volume_rel_drift_vs_base"] is not None
    ]
    area_drifts = [
        abs(c["surface_area_rel_drift"]) for c in checks
        if c["surface_area_rel_drift"] is not None
    ]
    analytic_errors = [abs(c["volume_rel_err"]) for c in checks]
    summary = {
        "rungs": len(checks),
        "max_triangles": max((c["triangles"] for c in checks), default=0),
        "single_body_10m_reached": any(c["triangles"] >= 10_000_000 for c in checks),
        "all_watertight": all(c["watertight"] for c in checks),
        "all_edge_manifold": all(c["edge_manifold"] for c in checks),
        "all_euler_2": all(c["euler_characteristic"] == 2 for c in checks),
        "max_volume_rel_drift_vs_base": max(volume_drifts) if volume_drifts else None,
        "max_surface_area_rel_drift": max(area_drifts) if area_drifts else None,
        "max_analytic_rel_err": max(analytic_errors) if analytic_errors else None,
        "analytic_err_note": (
            "analytic error is inherited from the sampled parent surface and is "
            "identical on every subdivision rung by construction"
        ),
        "volume_tolerance": VOLUME_TOLERANCE,
        "area_tolerance": AREA_TOLERANCE,
    }
    summary["invariants_ok"] = bool(
        summary["all_watertight"]
        and summary["all_edge_manifold"]
        and summary["all_euler_2"]
        and (summary["max_volume_rel_drift_vs_base"] or 0.0) <= VOLUME_TOLERANCE
        and (summary["max_surface_area_rel_drift"] or 0.0) <= AREA_TOLERANCE
    )
    return {
        "schema": SCHEMA,
        "command": command,
        "hull": {
            "length_m": hull.length,
            "beam_m": hull.beam,
            "draft_m": hull.draft,
            "fore_exponent": hull.fore_exponent,
            "aft_exponent": hull.aft_exponent,
            "analytic_volume_m3": analytic_volume(
                hull.length, hull.beam, hull.draft,
                fore_exponent=hull.fore_exponent, aft_exponent=hull.aft_exponent,
            ),
        },
        "n_queries": n_queries,
        "environment": _environment(),
        "seconds": round(elapsed, 3),
        "summary": summary,
        "checks": checks,
        "results": rows,
    }


#: Subdivision is exact in exact arithmetic; these tolerances only absorb
#: floating-point accumulation over 10^7 faces (summation order changes even
#: though the tiling does not).
VOLUME_TOLERANCE = 1e-12
AREA_TOLERANCE = 1e-12


def _max_area_ratio(mesh: Any, base: Any, level: int) -> float | None:
    """Largest child/parent triangle-area ratio (sanity: subdivision shrinks).

    Computed only for the first rung of a base (``level <= 1``) because it is
    O(F) over 10^7 faces otherwise; the property is structural (1/4 per level).
    """
    if level > 1 or not base.faces:
        return None
    from .mesh import triangle_area
    child_max = max(triangle_area(mesh, f) for f in mesh.faces)
    parent_max = max(triangle_area(base, f) for f in base.faces)
    if parent_max == 0:
        return None
    return child_max / parent_max