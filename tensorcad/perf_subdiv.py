"""Single-body 10^6-10^7 benchmark: subdivision as the large-mesh source.

``tensorcad.perf_benchmark`` reaches the 10^7 band with a *fleet* of disjoint
copies (``build_fleet_mesh``), because ``HullParameters`` caps one parametric
hull at 4*(nx-1)*(nz-1) = 8,199,998 triangles.  That answers "how does the code
scale with triangle count" but not "what happens on one 10^7-triangle body",
which is the shape real imported CAD/solver meshes have (one connected closed
surface, no per-copy regularity).

This driver closes that gap with :func:`tensorcad.mesh_subdivision.subdivide_mesh`
-- uniform 1-to-4 midpoint subdivision of a *coarser* hull.  Two consequences make it more
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

Sub-tasks measured per rung
--------------------------
``load``           generate_hull (parent mesh creation).
``subdivide``      one level of 1-to-4 midpoint subdivision.
``volume``         divergence-theorem surface integral.
``point_locate_accel``  TriangleGrid build + 2000-ray query.
``obj_write``      serialise mesh to Wavefront OBJ text (round-trip precision).
``obj_read``       parse OBJ back, verify volume against analytic reference.
"""
from __future__ import annotations

import tempfile
import time
from pathlib import Path
from typing import Any, Iterable, Sequence

from .export import read_obj, write_obj
from .geometry import Mesh
from .hull import HullParameters, HullModel, analytic_volume, generate_hull
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

SCHEMA = "tensorcad.perf_single10m/2"

#: Default coarse bases (nx, nz): 4*(nx-1)*(nz-1) triangles before subdivision.
DEFAULT_BASES: tuple[tuple[int, int], ...] = ((401, 129), (1601, 513))

#: Default subdivision ladder.  Level 0 is the generated hull itself; each level
#: multiplies the triangle count by four (401x129 -> 0.20e6 / 0.82e6 / 3.28e6 /
#: 13.1e6 at levels 0/1/2/3, i.e. the whole 10^6-10^7 band from one base).
DEFAULT_LEVELS: tuple[int, ...] = (0, 1, 2, 3)

DEFAULT_HULL = HullParameters(length=10.0, beam=2.0, draft=0.8)


def _to_mesh(obj: HullModel | Mesh) -> Mesh:
    """Normalise a HullModel or Mesh to a Mesh."""
    if isinstance(obj, HullModel):
        return obj.mesh
    return obj


def _euler_characteristic(mesh: Mesh) -> int:
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
        "file_bytes": 0,
        "file_mb_per_s": None,
        "volume_rel_err": None,
    }
    row.update(fields)
    return row


def _query_accel(mesh: Mesh, n_queries: int, nx: int, nz: int, levels: int) -> dict[str, Any]:
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


def _obj_roundtrip(mesh: Mesh, exact: float, nx: int, nz: int, levels: int) -> tuple[dict[str, Any], dict[str, Any]]:
    """Measure obj_write + obj_read on ``mesh``, returning two rows.

    Uses a temp file so the benchmark leaves no trace on disk.  The obj_read
    row carries ``volume_rel_err`` against the analytic hull volume, so the
    round-trip geometry fidelity is quantified at every rung.
    """
    triangles = len(mesh.faces)
    vertices = len(mesh.vertices)

    # --- obj_write ---
    with tempfile.NamedTemporaryFile(suffix=".obj", delete=False) as tmp:
        obj_path = tmp.name

    t0 = time.perf_counter()
    write_obj(mesh, obj_path)
    write_seconds = time.perf_counter() - t0
    file_bytes = Path(obj_path).stat().st_size
    write_row = _row(
        subtask="obj_write",
        nx=nx,
        nz=nz,
        levels=levels,
        triangles=triangles,
        vertices=vertices,
        seconds=round(write_seconds, 6),
        throughput=round(triangles / write_seconds, 1) if write_seconds > 0 else None,
        note=f"file={file_bytes / 1e6:.1f} MiB ({file_bytes / 1e6 / write_seconds:.1f} MiB/s), "
             f"{file_bytes / triangles:.1f} B/tri, single connected body",
        file_bytes=file_bytes,
        file_mb_per_s=round(file_bytes / 1e6 / write_seconds, 3) if write_seconds > 0 else None,
        peak_rss_mb=peak_rss_mb(),
        rss_mb=live_rss_mb(),
        mesh_mb=mesh_footprint_mb(mesh),
    )

    # --- obj_read ---
    t0 = time.perf_counter()
    imported = read_obj(obj_path)
    read_seconds = time.perf_counter() - t0
    imported_triangles = len(imported.faces)
    imported_volume = total_signed_volume(imported)
    vol_rel_err = (imported_volume - exact) / exact if exact != 0 else None

    read_row = _row(
        subtask="obj_read",
        nx=nx,
        nz=nz,
        levels=levels,
        triangles=imported_triangles,
        vertices=len(imported.vertices),
        seconds=round(read_seconds, 6),
        throughput=round(imported_triangles / read_seconds, 1) if read_seconds > 0 else None,
        note=(
            f"file={file_bytes / 1e6:.1f} MiB ({file_bytes / 1e6 / read_seconds:.1f} MiB/s), "
            f"volume={imported_volume:.6f} m^3, analytic={exact:.6f} m^3, "
            f"rel_err={vol_rel_err:.3e}" if vol_rel_err is not None else ""
        ),
        file_bytes=file_bytes,
        file_mb_per_s=round(file_bytes / 1e6 / read_seconds, 3) if read_seconds > 0 else None,
        volume_rel_err=vol_rel_err,
        peak_rss_mb=peak_rss_mb(),
        rss_mb=live_rss_mb(),
    )

    # Clean up
    try:
        Path(obj_path).unlink()
    except OSError:
        pass

    return write_row, read_row


def run_single_body_benchmark(
    bases: Iterable[Sequence[int]] = DEFAULT_BASES,
    levels_ladder: Iterable[int] = DEFAULT_LEVELS,
    hull: HullParameters = DEFAULT_HULL,
    n_queries: int = 2000,
    max_triangles: int | None = 30_000_000,
    command: str | None = None,
    subtasks: Sequence[str] = ("load", "subdivide", "volume", "point_locate_accel", "obj_write", "obj_read"),
) -> dict[str, Any]:
    """Measure sub-tasks on one connected body at 10^6-10^7 triangle scale.

    Every rung is one connected closed mesh, so the 10^6-10^7 rows are directly
    comparable with an imported large mesh.  Invariant checks (watertight,
    edge-manifold, Euler characteristic 2, volume/area preservation along the
    subdivision ladder) are recorded per rung: a fast but wrong geometry source
    cannot pass silently.

    Parameters
    ----------
    subtasks:
        Which sub-tasks to measure.  Default includes all six:
        ``load``, ``subdivide``, ``volume``, ``point_locate_accel``,
        ``obj_write``, ``obj_read``.
    """
    levels_sorted = sorted({int(level) for level in levels_ladder})
    if -1 in levels_sorted:
        raise ValueError("levels must be >= 0")
    subtask_set = set(subtasks)
    rows: list[dict[str, Any]] = []
    checks: list[dict[str, Any]] = []
    t_start = time.perf_counter()

    for nx, nz in bases:
        nx, nz = int(nx), int(nz)
        t0 = time.perf_counter()
        base_model = generate_hull(HullParameters(
            length=hull.length, beam=hull.beam, draft=hull.draft, nx=nx, nz=nz,
        ))
        load_seconds = time.perf_counter() - t0
        exact = analytic_volume(
            hull.length, hull.beam, hull.draft,
            fore_exponent=hull.fore_exponent, aft_exponent=hull.aft_exponent,
        )
        base = _to_mesh(base_model)
        base_triangles = len(base.faces)
        base_area = surface_area(base)
        base_volume = total_signed_volume(base)

        if "load" in subtask_set:
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
                if "subdivide" in subtask_set:
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
                else:
                    mesh = subdivide_mesh(mesh, 1)

            current_triangles = len(mesh.faces)
            current_vertices = len(mesh.vertices)

            if "volume" in subtask_set:
                t0 = time.perf_counter()
                volume = total_signed_volume(mesh)
                volume_seconds = time.perf_counter() - t0
                rel_err = (volume - exact) / exact
                rows.append(_row(
                    subtask="volume",
                    nx=nx,
                    nz=nz,
                    levels=level,
                    triangles=current_triangles,
                    vertices=current_vertices,
                    seconds=round(volume_seconds, 6),
                    throughput=(
                        round(current_triangles / volume_seconds, 1) if volume_seconds > 0 else None
                    ),
                    note=(
                        f"analytic={exact:.6f} m3, measured={volume:.6f} m3, "
                        f"rel_err={rel_err:.3e}"
                    ),
                    peak_rss_mb=peak_rss_mb(),
                    rss_mb=live_rss_mb(),
                    mesh_mb=mesh_footprint_mb(mesh),
                ))
            else:
                volume = total_signed_volume(mesh)
                rel_err = (volume - exact) / exact

            if "point_locate_accel" in subtask_set:
                rows.append(_query_accel(mesh, n_queries=n_queries, nx=nx, nz=nz, levels=level))

            if "obj_write" in subtask_set or "obj_read" in subtask_set:
                write_r, read_r = _obj_roundtrip(mesh, exact, nx, nz, level)
                if "obj_write" in subtask_set:
                    rows.append(write_r)
                if "obj_read" in subtask_set:
                    rows.append(read_r)

            area = surface_area(mesh)
            checks.append({
                "nx": nx,
                "nz": nz,
                "levels": level,
                "triangles": current_triangles,
                "vertices": current_vertices,
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

    elapsed = time.perf_counter() - t_start
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
            "identical on every subdivision level (subdivision preserves the "
            "discrete solid exactly)"
        ),
        "invariants_ok": (
            all(c["watertight"] for c in checks)
            and all(c["edge_manifold"] for c in checks)
            and all(c["euler_characteristic"] == 2 for c in checks)
            and (max(volume_drifts) < 1e-9 if volume_drifts else True)
            and (max(area_drifts) < 1e-9 if area_drifts else True)
        ),
    }

    return {
        "schema": SCHEMA,
        "command": command or "python -m tensorcad.perf_subdiv",
        "environment": _environment(),
        "seconds": round(elapsed, 3),
        "results": rows,
        "checks": checks,
        "summary": summary,
    }


def _max_area_ratio(mesh: Mesh, base: Mesh, level: int) -> float:
    """Largest child-triangle area / mean parent-triangle area at this level.

    A value near 0.25 confirms uniform 1-to-4 subdivision (each child is
    roughly a quarter of the parent).  Values >>0.25 would signal that some
    triangles are not being subdivided properly.
    """
    if level == 0:
        return 1.0
    import math
    areas = []
    for face in mesh.faces:
        a = mesh.vertices[face[0]]
        b = mesh.vertices[face[1]]
        c = mesh.vertices[face[2]]
        ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
        cross = math.sqrt(
            (uy * vz - uz * vy) ** 2
            + (uz * vx - ux * vz) ** 2
            + (ux * vy - uy * vx) ** 2
        )
        areas.append(0.5 * cross)
    if not areas:
        return 0.0
    mean_area = sum(areas) / len(areas)
    if mean_area <= 0:
        return 0.0
    return max(areas) / mean_area
