"""Large-scale geometry performance benchmark (10^5--10^7 triangle regime).

Answers the question in GAP-CAD-D4-PERF-SCALE: *at what mesh size does each
geometry operation become expensive, and when does it need a cluster?*

Four sub-tasks are measured on the same mesh ladder:

``load``          Mesh construction from analytic parameters (hull generation).
``point_locate``  Ray-casting / signed-distance queries, brute-force reference.
``point_locate_accel``  Same queries through :class:`tensorcad.spatial_index.TriangleGrid`
                  (uniform grid + DDA), the path that makes 10^6--10^7 practical.
``voxelize``      SDF sampling onto a regular grid (to_voxel_grid).
``volume``        Surface-integral volume by the divergence theorem.
``obj_write``     Serialising the rung's mesh to Wavefront OBJ text (file export).
``obj_read``      Parsing that OBJ file back into a mesh (file import), with an
                  as-imported volume check against the analytic reference.

Each sub-task is run at several mesh sizes (triangle counts).  The benchmark
records wall-clock seconds, triangle count, and a derived throughput metric
(triangles/s or queries/s).  A *regression gate* compares each run against a
stored baseline; exceeding 2x the baseline time flags RED.

The module is standard-library only so it can run in any TensorCAD environment.
NumPy is imported lazily inside the voxel path.

Usage
-----
::

    # Quick smoke (small meshes, ~1 s total)
    python -m tensorcad.perf_benchmark --quick

    # Full ladder (up to ~10^7 triangles, may take minutes)
    python -m tensorcad.perf_benchmark --full

    # Save baseline for future regression gates (write it under benchmarks/,
    # which is version-controlled; ``runs/`` is filtered out by the merge gate,
    # so a baseline saved there never reaches ``main`` and the gate goes blind)
    python -m tensorcad.perf_benchmark --quick --save-baseline benchmarks/perf_scale/baseline.json

    # Check against the committed baseline
    python -m tensorcad.perf_benchmark --quick --baseline benchmarks/perf_scale/baseline.json

See ``docs/perf_benchmark.md`` for the gate conventions (this module's
documentation of record) and ``docs/perf_scale_10m.md`` for the 10^6--10^7
measurements.

Exit code is 0 when all gates pass, 1 when any gate fails.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import shutil
import sys
import tempfile
import time
from collections import OrderedDict
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Data contracts
# ---------------------------------------------------------------------------

BENCHMARK_SCHEMA = "tensorcad.perf_benchmark/2"

#: Status tags for individual measurements.
MEASURED = "measured"
NOT_RUN = "not_run"
FAILED = "failed"

#: The mesh ladder: (nx, nz) pairs that produce increasing triangle counts.
#: The default ladder goes from ~2k to ~10M triangles.
#: Triangle count of a rung is ``4 * (nx - 1) * (nz - 1)`` (two 4-sided strips
#: per grid cell); the comments below are the measured/estimated counts, not
#: guesses, and the top rung lands at 1.0e7 triangles.
DEFAULT_LADDER: tuple[tuple[int, int], ...] = (
    (21, 9),      # 640 triangles
    (41, 13),     # 1,920
    (81, 21),     # 6,400
    (161, 41),    # 25,600
    (321, 81),    # 102,400
    (641, 161),   # 409,600
    (1281, 321),  # 1,638,400  -> 10^6 band
    (2001, 501),  # 4,000,000  -> mid 10^6 band
    (4001, 513),  # 8,192,000 -> top of the constructible ladder; HullParameters
                  # caps nx at 4001 stations and nz at 513 levels, so 8.19e6 is
                  # the largest parametric hull (a 10^7 rung needs an imported
                  # mesh rather than generate_hull)
)

#: Quick ladder for CI / smoke tests.
QUICK_LADDER: tuple[tuple[int, int], ...] = (
    (21, 9),
    (41, 13),
    (81, 21),
    (161, 41),
)

#: Voxel grid resolution ladder (scales with mesh complexity).
VOXEL_LADDER: tuple[int, ...] = (32, 48, 64, 96, 128, 192, 256, 384)
QUICK_VOXEL_LADDER: tuple[int, ...] = (32, 48, 64, 96)

#: Maximum triangles for a single run before we warn about cluster need.
#: This is a *warning* threshold for the naive path, not a measured cluster
#: boundary: the accelerated path runs 12e6 triangles in seconds, so the real
#: cluster question is algorithmic (see ``docs/perf_scale_10m.md`` §5).
CLUSTER_THRESHOLD = 5_000_000

#: Default RED threshold for the wall-clock regression gate: a measurement that
#: takes more than this multiple of the stored baseline seconds is a regression.
DEFAULT_GATE_FACTOR = 2.0

#: Ray queries for point_locate (scales with mesh size in full mode).
QUICK_N_QUERIES = 1000
FULL_N_QUERIES = 10000


@dataclass(frozen=True)
class SubTaskResult:
    """One measurement of one sub-task at one mesh/grid size."""

    subtask: str          # "load" | "point_locate" | "voxelize" | "volume"
    nx: int               # For mesh tasks: nx; for voxel: grid resolution
    nz: int               # For mesh tasks: nz; for voxel: 0
    triangles: int        # Mesh triangle count (0 for voxel-only)
    vertices: int         # Mesh vertex count (0 for voxel-only)
    grid_cells: int       # Voxel grid cells (0 for mesh-only)
    seconds: float | None
    status: str = MEASURED
    throughput: float | None = None   # triangles/s, queries/s, or cells/s
    note: str = ""
    #: Index build time for accelerated sub-tasks (spatial index construction).
    build_seconds: float | None = None
    #: Spatial-index size descriptors (0 for non-indexed sub-tasks).
    index_cells: int = 0
    index_entries: int = 0
    index_bytes: int = 0
    #: Peak resident set size of the process when this measurement finished.
    #: This is ``ru_maxrss``, a *process-wide high-water mark* that cannot be
    #: reset between rungs: read it as "process peak so far", never as "this
    #: rung's peak".  :attr:`rss_mb` is the live footprint at rung end and does
    #: distinguish rungs.
    peak_rss_mb: float | None = None
    #: Live resident set size (``/proc/self/status`` ``VmRSS``) at rung end, MiB.
    rss_mb: float | None = None
    #: Estimated bytes->MiB footprint of the measured mesh itself (sampled
    #: ``sys.getsizeof`` extrapolation, see :func:`mesh_footprint_mb`).  This is
    #: the memory question the scale table has to answer ("how much RAM do 10^7
    #: triangles need"), and unlike RSS it does not depend on machine load.
    mesh_mb: float | None = None
    #: Number of disjoint bodies in the measured mesh.  ``1`` is a single
    #: parametric hull; ``>1`` marks a multi-body *fleet* rung, which exists
    #: because ``HullParameters`` caps one body at 8.19e6 triangles and the
    #: 10^7 band is otherwise unreachable (see :func:`build_fleet_mesh`).
    copies: int = 1
    #: Size of the OBJ file the import sub-tasks wrote/read (0 elsewhere).  A
    #: structured field rather than text in ``note`` so the published table can
    #: be rendered from the artifact instead of re-typing numbers.
    file_bytes: int = 0
    #: File throughput of ``obj_write``/``obj_read`` (MiB/s).
    file_mb_per_s: float | None = None
    #: Signed-volume error of the *as-imported* mesh against the analytic
    #: reference (``obj_read``); the parse is only a measurement if the geometry
    #: that came back is still the geometry that went in.
    volume_rel_err: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "subtask": self.subtask,
            "copies": self.copies,
            "nx": self.nx,
            "nz": self.nz,
            "triangles": self.triangles,
            "vertices": self.vertices,
            "grid_cells": self.grid_cells,
            "seconds": self.seconds,
            "status": self.status,
            "throughput": self.throughput,
            "note": self.note,
            "build_seconds": self.build_seconds,
            "index_cells": self.index_cells,
            "index_entries": self.index_entries,
            "index_bytes": self.index_bytes,
            "peak_rss_mb": self.peak_rss_mb,
            "rss_mb": self.rss_mb,
            "mesh_mb": self.mesh_mb,
            "file_bytes": self.file_bytes,
            "file_mb_per_s": self.file_mb_per_s,
            "volume_rel_err": self.volume_rel_err,
        }


@dataclass(frozen=True)
class GateVerdict:
    """Pass/fail for one regression gate."""

    name: str
    passed: bool
    detail: str = ""
    baseline_seconds: float | None = None
    measured_seconds: float | None = None
    ratio: float | None = None
    #: ``measured`` when the gate had enough data to judge, ``not_run`` when it
    #: did not.  A ``not_run`` gate is not a pass: it is an explicit statement
    #: that the check could not be performed on the data at hand.
    status: str = MEASURED

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "passed": self.passed,
            "status": self.status,
            "detail": self.detail,
            "baseline_seconds": self.baseline_seconds,
            "measured_seconds": self.measured_seconds,
            "ratio": self.ratio,
        }


@dataclass(frozen=True)
class PerfBenchmarkResult:
    """Complete performance benchmark result."""

    results: tuple[SubTaskResult, ...]
    gates: tuple[GateVerdict, ...] = ()
    environment: dict[str, Any] = field(default_factory=dict)
    seconds: float = 0.0
    command: str = ""
    cluster_warnings: tuple[str, ...] = ()

    @property
    def passed(self) -> bool:
        return all(g.passed for g in self.gates)

    @property
    def gates_checked(self) -> int:
        """Number of gates that actually produced a verdict.

        A ``not_run`` gate carries ``passed=True`` for backwards compatibility,
        so ``passed`` alone cannot tell "everything was checked and green" from
        "nothing was checked".  This count can (see :attr:`verdict`).
        """
        return sum(1 for g in self.gates if g.status != NOT_RUN)

    @property
    def gates_not_run(self) -> int:
        return sum(1 for g in self.gates if g.status == NOT_RUN)

    @property
    def verdict(self) -> str:
        """``fail`` / ``pass`` / ``not_run`` - honest three-way verdict.

        ``not_run`` means no gate could be evaluated (e.g. a ladder with fewer
        than ``MIN_SCALING_POINTS`` rungs): that is *not* a pass.
        """
        if any(not g.passed for g in self.gates):
            return "fail"
        return "pass" if self.gates_checked else "not_run"

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": BENCHMARK_SCHEMA,
            "passed": self.passed,
            "verdict": self.verdict,
            "gates_checked": self.gates_checked,
            "gates_not_run": self.gates_not_run,
            "seconds": self.seconds,
            "command": self.command,
            "environment": self.environment,
            "cluster_warnings": list(self.cluster_warnings),
            "gates": [g.to_dict() for g in self.gates],
            "results": [r.to_dict() for r in self.results],
        }

    def write_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_dict(), indent=2) + "\n", encoding="utf-8")
        return target


# ---------------------------------------------------------------------------
# Environment capture
# ---------------------------------------------------------------------------

def _environment() -> dict[str, Any]:
    from tensorcad import __version__
    return {
        "tensorcad": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "loadavg": _loadavg(),
    }


def _loadavg() -> dict[str, float]:
    try:
        parts = open("/proc/loadavg").read().split()
        return {"load1": float(parts[0]), "load5": float(parts[1]), "load15": float(parts[2])}
    except Exception:
        return {}


# ---------------------------------------------------------------------------
# Shared measurement helpers
# ---------------------------------------------------------------------------

def peak_rss_mb() -> float | None:
    """Peak resident set size of this process in MiB (None if unavailable).

    Recorded per measurement so the scale table answers "how much memory does
    10^7 triangles need", not only "how many seconds".
    """
    try:
        import resource

        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    except Exception:  # pragma: no cover - non-POSIX
        return None
    # Linux reports KiB, macOS reports bytes.
    if sys.platform == "darwin":  # pragma: no cover
        return rss / (1024.0 * 1024.0)
    return rss / 1024.0


def live_rss_mb() -> float | None:
    """Current resident set size in MiB (Linux ``VmRSS``), ``None`` elsewhere."""
    try:
        with open("/proc/self/status", encoding="ascii", errors="replace") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024.0
    except Exception:  # pragma: no cover - non-Linux
        return None
    return None  # pragma: no cover


def mesh_footprint_mb(mesh: Any, sample: int = 256) -> float | None:
    """Estimate the memory held by a triangle mesh, in MiB.

    ``sys.getsizeof`` over every vertex of a 10^7-triangle mesh would itself
    cost seconds, so the per-element cost is measured on at most ``sample``
    evenly spaced elements and extrapolated.  Vertices are tuples of floats and
    faces are tuples of flat vertex indices (``tensorcad.geometry.Mesh``), so
    the two are scored with different element sizes.

    The result is *mesh-only* (the geometry containers, not the spatial index
    and not the interpreter), which is what "how much RAM do 10^7 triangles
    need" actually asks.  It is an estimate by construction: small ints are
    CPython-interned in the low range, so face indices are charged the
    measured sample's average, not a per-index worst case.
    """
    try:
        import sys as _sys

        verts = getattr(mesh, "vertices", None)
        faces = getattr(mesh, "faces", None)
        if verts is None or faces is None:
            return None
        nv, nf = len(verts), len(faces)
        if nv == 0 or nf == 0:
            return None

        vs = max(1, min(sample, nv))
        vstep = max(1, nv // vs)
        vsum, vcnt = 0, 0
        for i in range(0, nv, vstep):
            vertex = verts[i]
            vsum += _sys.getsizeof(vertex) + sum(_sys.getsizeof(c) for c in vertex)
            vcnt += 1
            if vcnt >= vs:
                break

        fs = max(1, min(sample, nf))
        fstep = max(1, nf // fs)
        fsum, fcnt = 0, 0
        for i in range(0, nf, fstep):
            face = faces[i]
            fsum += _sys.getsizeof(face) + sum(_sys.getsizeof(k) for k in face)
            fcnt += 1
            if fcnt >= fs:
                break

        per_v = vsum / vcnt
        per_f = fsum / fcnt
        total = nv * per_v + nf * per_f + _sys.getsizeof(verts) + _sys.getsizeof(faces)
        return total / (1024.0 * 1024.0)
    except Exception:  # pragma: no cover - defensive
        return None


#: Bound on retained meshes: one rung plus its predecessor.
MESH_CACHE_SIZE = 2
_MESH_CACHE: "OrderedDict[tuple[int, int], Any]" = OrderedDict()


def _put_mesh(nx: int, nz: int, model: Any) -> None:
    _MESH_CACHE[(nx, nz)] = model
    _MESH_CACHE.move_to_end((nx, nz))
    while len(_MESH_CACHE) > MESH_CACHE_SIZE:
        _MESH_CACHE.popitem(last=False)


def _cached_hull(nx: int, nz: int):
    """Generate (and keep) one hull mesh per resolution.

    The ladder regenerates the same mesh for every sub-task, which makes the
    10^6--10^7 rungs spend all their wall clock in the generator and is simply
    unaffordable once a single rung costs minutes.
    """
    key = (nx, nz)
    model = _MESH_CACHE.get(key)
    if model is None:
        from tensorcad.hull import HullParameters, generate_hull

        model = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=nx, nz=nz))
        _put_mesh(nx, nz, model)
    else:
        _MESH_CACHE.move_to_end(key)
    return model


def clear_mesh_cache() -> None:
    """Drop cached meshes (tests and repeated benchmark runs)."""
    _MESH_CACHE.clear()

# ---------------------------------------------------------------------------
# Multi-body source for the 10^7 band
# ---------------------------------------------------------------------------

#: Centre-to-centre spacing between hull copies, in beams.  One beam of gap
#: keeps the bodies strictly disjoint, so the assembled mesh has no shared
#: faces and the signed-volume integral stays exact (``copies * V``).
FLEET_SPACING_BEAMS = 2.0

#: Bound on retained fleet meshes: one rung at a time.  A 10^7-triangle fleet
#: is several GB of Python tuples, keeping two would OOM the node.
_FLEET_CACHE: "OrderedDict[tuple[int, int, int], Any]" = OrderedDict()
FLEET_CACHE_SIZE = 1


def build_fleet_mesh(model: Any, copies: int, spacing: float | None = None) -> Any:
    """Assemble ``copies`` translated copies of one hull into a single mesh.

    Why this exists: ``HullParameters`` refuses ``nx > 4001`` stations and
    ``nz > 513`` levels, so the largest *parametric* hull is 8,192,000
    triangles and the 10^7 band of GAP-CAD-D4-PERF-SCALE cannot be reached by
    ``generate_hull`` alone.  Translating k disjoint copies gives a k-fold
    larger mesh at O(N) copy cost (no generator loop), and because the bodies
    do not touch, every geometry operation used here keeps a closed-form
    answer: total signed volume is exactly ``copies`` x the single hull volume.

    The mesh is *not* a substitute for a single large imported body when the
    question is "how does the generator scale" - the assembly cost is measured
    and reported separately as the ``fleet_assemble`` sub-task for exactly that
    reason.  The brute-force ``_RayCaster`` reference is star-shaped only and
    must not be used on a fleet.

    Parameters
    ----------
    model:
        A generated hull model (``model.mesh`` is the :class:`~tensorcad.geometry.Mesh`).
    copies:
        Number of bodies, ``>= 1``.
    spacing:
        Centre-to-centre spacing in model units; defaults to
        ``FLEET_SPACING_BEAMS * beam`` of the hull parameters.
    """
    if copies < 1:
        raise ValueError("copies must be >= 1")

    from tensorcad.geometry import Mesh

    base = model.mesh
    if copies == 1:
        return base
    if spacing is None:
        beam = getattr(getattr(model, "parameters", None), "beam", None)
        if beam is None:
            raise ValueError("spacing is required when the model has no beam parameter")
        spacing = FLEET_SPACING_BEAMS * float(beam)

    base_verts = base.vertices
    base_faces = base.faces
    n_verts = len(base_verts)
    vertices: list[Any] = []
    faces: list[Any] = []
    for k in range(copies):
        dy = k * spacing
        vertices.extend((x, y + dy, z) for (x, y, z) in base_verts)
        offset = k * n_verts
        faces.extend((a + offset, b + offset, c + offset) for (a, b, c) in base_faces)
    return Mesh(vertices=tuple(vertices), faces=tuple(faces), units=base.units,
                coordinate_frame=base.coordinate_frame)


def _cached_fleet(nx: int, nz: int, copies: int) -> Any:
    """Return (and cache) the fleet mesh for one (resolution, body count) key."""
    key = (nx, nz, copies)
    mesh = _FLEET_CACHE.get(key)
    if mesh is None:
        mesh = build_fleet_mesh(_cached_hull(nx, nz), copies)
        _FLEET_CACHE[key] = mesh
        _FLEET_CACHE.move_to_end(key)
        while len(_FLEET_CACHE) > FLEET_CACHE_SIZE:
            _FLEET_CACHE.popitem(last=False)
    else:
        _FLEET_CACHE.move_to_end(key)
    return mesh


def clear_fleet_cache() -> None:
    """Drop cached fleet meshes."""
    _FLEET_CACHE.clear()


def _mesh_for(nx: int, nz: int, copies: int) -> Any:
    """Mesh accessor shared by the mesh sub-tasks (single body or fleet)."""
    return _cached_hull(nx, nz).mesh if copies == 1 else _cached_fleet(nx, nz, copies)


# ---------------------------------------------------------------------------
# Sub-task implementations
# ---------------------------------------------------------------------------

def _run_load(nx: int, nz: int) -> SubTaskResult:
    """Measure hull generation time (always a fresh generation).

    The generated mesh is handed to the shared cache so the following
    sub-tasks at the same rung do not regenerate it; the timing here covers
    generation only.
    """
    from tensorcad.hull import HullParameters, generate_hull
    params = HullParameters(length=10.0, beam=2.0, draft=0.8, nx=nx, nz=nz)
    t0 = time.perf_counter()
    model = generate_hull(params)
    elapsed = time.perf_counter() - t0
    _put_mesh(nx, nz, model)
    tris = len(model.mesh.faces)
    verts = len(model.mesh.vertices)
    return SubTaskResult(
        subtask="load", nx=nx, nz=nz, triangles=tris, vertices=verts,
        grid_cells=0,
        seconds=round(elapsed, 6),
        throughput=round(tris / elapsed, 1) if elapsed > 0 else None,
        peak_rss_mb=peak_rss_mb(),
        rss_mb=live_rss_mb(),
        mesh_mb=mesh_footprint_mb(model.mesh),
    )


def _run_point_locate(nx: int, nz: int, n_queries: int = 1000, copies: int = 1) -> SubTaskResult:
    """Measure ray-casting / point-location throughput.

    The _RayCaster casts rays from the origin; the hull is centered at the
    origin so rays in all directions that intersect the hull will return a
    radius.  We use a Fibonacci sphere for uniform directional coverage.
    """
    from tensorcad.benchmark import _RayCaster

    if copies != 1:
        raise ValueError(
            "_RayCaster is a star-shaped single-body reference; "
            "use point_locate_accel for multi-body fleet rungs"
        )
    model = _cached_hull(nx, nz)
    tris = len(model.mesh.faces)
    verts = len(model.mesh.vertices)

    directions = _fibonacci_directions(n_queries)
    caster = _RayCaster(model.mesh)

    t0 = time.perf_counter()
    hits = 0
    for direction in directions:
        hit = caster.radius(direction)
        if hit is not None:
            hits += 1
    elapsed = time.perf_counter() - t0

    return SubTaskResult(
        subtask="point_locate", nx=nx, nz=nz, triangles=tris, vertices=verts,
        grid_cells=0,
        seconds=round(elapsed, 6),
        throughput=round(n_queries / elapsed, 1) if elapsed > 0 else None,
        note=f"hits={hits}/{n_queries}, brute-force reference (no index)",
        peak_rss_mb=peak_rss_mb(),
        rss_mb=live_rss_mb(),
        mesh_mb=mesh_footprint_mb(model.mesh),
    )


def _run_point_locate_accel(
    nx: int, nz: int, n_queries: int = 1000, copies: int = 1
) -> SubTaskResult:
    """Measure the indexed point-location path (uniform grid + DDA).

    ``seconds`` is the query time only; ``build_seconds`` is reported
    separately because index construction is a one-off cost.  The two numbers
    answer different engineering questions ("can I afford one index build" vs
    "what does a query cost").

    With ``copies > 1`` the measured mesh is a multi-body fleet (see
    :func:`build_fleet_mesh`); the queries are the same Fibonacci directions
    from the origin, so for a fleet they report the nearest body's intersection.
    """
    from tensorcad.spatial_index import TriangleGrid

    mesh = _mesh_for(nx, nz, copies)
    tris = len(mesh.faces)
    verts = len(mesh.vertices)

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

    return SubTaskResult(
        subtask="point_locate_accel", nx=nx, nz=nz, triangles=tris, vertices=verts,
        copies=copies,
        grid_cells=int(stats["cells"]),
        seconds=round(elapsed, 6),
        throughput=round(n_queries / elapsed, 1) if elapsed > 0 else None,
        build_seconds=round(build, 6),
        index_cells=int(stats["cells"]),
        index_entries=int(stats["entries"]),
        index_bytes=int(stats["index_bytes"]),
        peak_rss_mb=peak_rss_mb(),
        rss_mb=live_rss_mb(),
        mesh_mb=mesh_footprint_mb(mesh),
        note=(
            f"hits={hits}/{n_queries}, grid={tuple(stats['resolution'])}, "
            f"entries/tri={stats['entries_per_triangle']:.2f}, "
            f"build={build:.3f}s, bodies={copies}"
        ),
    )


def _run_voxelize(grid_res: int) -> SubTaskResult:
    """Measure SDF voxelization time using SuboffHullField at given grid resolution."""
    try:
        import numpy as np  # noqa: F401
    except ImportError:
        return SubTaskResult(
            subtask="voxelize", nx=grid_res, nz=0, triangles=0, vertices=0,
            grid_cells=0, seconds=None, status=NOT_RUN, note="numpy not available",
        )

    from tensorcad.suboff import SuboffHullField
    from tensorcad.voxel import to_voxel_grid

    field = SuboffHullField()

    t0 = time.perf_counter()
    grid = to_voxel_grid(field, grid_res)
    elapsed = time.perf_counter() - t0

    cells = grid.cells
    return SubTaskResult(
        subtask="voxelize", nx=grid_res, nz=0, triangles=0, vertices=0,
        grid_cells=cells,
        seconds=round(elapsed, 6),
        throughput=round(cells / elapsed, 1) if elapsed > 0 else None,
        note=f"grid={grid.shape}, cells={cells:,}",
        peak_rss_mb=peak_rss_mb(),
        rss_mb=live_rss_mb(),
    )


def _run_volume(nx: int, nz: int, copies: int = 1) -> SubTaskResult:
    """Measure volume integration time (divergence theorem over the triangles).

    The measured volume is reported together with its analytic reference
    (``copies`` x ``hull.analytic_volume``) and the relative error, so a rung
    that produces a fast-but-wrong number cannot pass unnoticed.
    """
    from tensorcad import mesh as meshlib
    from tensorcad.hull import analytic_volume

    mesh = _mesh_for(nx, nz, copies)
    tris = len(mesh.faces)
    verts = len(mesh.vertices)

    t0 = time.perf_counter()
    vol = meshlib.total_signed_volume(mesh)
    elapsed = time.perf_counter() - t0

    reference = copies * analytic_volume(10.0, 2.0, 0.8)
    rel_err = (vol - reference) / reference if reference else float("nan")

    return SubTaskResult(
        subtask="volume", nx=nx, nz=nz, triangles=tris, vertices=verts,
        grid_cells=0,
        copies=copies,
        seconds=round(elapsed, 6),
        throughput=round(tris / elapsed, 1) if elapsed > 0 else None,
        note=(
            f"volume={vol:.6f} m^3, analytic={reference:.6f} m^3, "
            f"rel_err={rel_err:+.3%}, bodies={copies}"
        ),
        peak_rss_mb=peak_rss_mb(),
        rss_mb=live_rss_mb(),
        mesh_mb=mesh_footprint_mb(mesh),
    )


def _run_fleet_assemble(nx: int, nz: int, copies: int) -> SubTaskResult:
    """Measure the O(N) multi-body assembly of a fleet mesh.

    Timed separately from ``load`` (generation) because they answer different
    questions: generation is a generator-loop cost, assembly is a pure copy.  A
    fleet rung is only interpretable if both numbers are on the table.
    """
    from tensorcad.mesh import total_signed_volume  # noqa: F401  (import check)

    base = _cached_hull(nx, nz)
    base_tris = len(base.mesh.faces)

    t0 = time.perf_counter()
    mesh = build_fleet_mesh(base, copies)
    elapsed = time.perf_counter() - t0

    _FLEET_CACHE[(nx, nz, copies)] = mesh
    _FLEET_CACHE.move_to_end((nx, nz, copies))
    while len(_FLEET_CACHE) > FLEET_CACHE_SIZE:
        _FLEET_CACHE.popitem(last=False)

    tris = len(mesh.faces)
    return SubTaskResult(
        subtask="fleet_assemble", nx=nx, nz=nz, triangles=tris, vertices=len(mesh.vertices),
        grid_cells=0,
        copies=copies,
        seconds=round(elapsed, 6),
        throughput=round(tris / elapsed, 1) if elapsed > 0 else None,
        note=f"bodies={copies} x {base_tris:,} tri from one generated body",
        peak_rss_mb=peak_rss_mb(),
        rss_mb=live_rss_mb(),
        mesh_mb=mesh_footprint_mb(mesh),
    )


# ---------------------------------------------------------------------------
# File import / export sub-tasks (OBJ at scale)
# ---------------------------------------------------------------------------

#: Directory the import sub-tasks write their OBJ file to.  Under the system
#: temp dir on purpose: one rung at 10^7 triangles is a ~1 GiB text file and it
#: must never land inside the repository.
IMPORT_DIR = Path(tempfile.gettempdir()) / "tensorcad-perf-import"

#: Bytes-per-triangle estimate of the OBJ layout written by ``write_obj`` (one
#: ``v`` line per vertex, one ``f`` line per face, ~1:1 vertices:faces on the
#: hull).  Only used for the free-space guard; the measured file size is always
#: what gets reported.
OBJ_BYTES_PER_TRIANGLE = 120

#: Free space is demanded as a multiple of the estimate so a write cannot fill
#: the filesystem up to the last byte.
OBJ_FREE_SPACE_FACTOR = 4.0


def import_obj_path(nx: int, nz: int, copies: int = 1) -> Path:
    """Path of the OBJ file the import sub-tasks share at one rung."""
    return IMPORT_DIR / f"hull_{nx}x{nz}_b{copies}.obj"


def clear_import_files() -> list[str]:
    """Delete every OBJ file written by the import sub-tasks."""
    removed: list[str] = []
    if IMPORT_DIR.is_dir():
        for path in sorted(IMPORT_DIR.glob("hull_*x*_b*.obj")):
            removed.append(path.name)
            path.unlink(missing_ok=True)
    return removed


def _prune_import_files(keep: Path) -> None:
    """Keep only ``keep`` so a ladder does not leave every rung on disk."""
    if not IMPORT_DIR.is_dir():
        return
    for path in IMPORT_DIR.glob("hull_*x*_b*.obj"):
        if path != keep:
            path.unlink(missing_ok=True)


def _free_mb(path: Path) -> float:
    probe = path
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    return shutil.disk_usage(probe).free / (1024 * 1024)


def _obj_space_guard(tris: int) -> str | None:
    """Reason string when there is not enough free space to write the file."""
    needed_mb = tris * OBJ_BYTES_PER_TRIANGLE / (1024 * 1024)
    free_mb = _free_mb(IMPORT_DIR)
    if free_mb < OBJ_FREE_SPACE_FACTOR * needed_mb:
        return (
            f"not enough free space under {IMPORT_DIR}: {free_mb:,.0f} MiB free, "
            f"~{needed_mb:,.0f} MiB needed for {tris:,} triangles"
        )
    return None


def _run_obj_write(nx: int, nz: int, copies: int = 1) -> SubTaskResult:
    """Measure writing the rung's mesh as OBJ text (serialisation cost).

    Deliberately split from :func:`_run_obj_read`: "importing this hull takes
    eight minutes" is not actionable, "writing costs X, parsing costs Y" is.
    The file is left on disk so the ``obj_read`` rung at the same size measures
    parsing only.
    """
    from tensorcad.export import write_obj

    mesh = _mesh_for(nx, nz, copies)
    tris = len(mesh.faces)
    path = import_obj_path(nx, nz, copies)
    reason = _obj_space_guard(tris)
    if reason is not None:
        return SubTaskResult(
            subtask="obj_write", nx=nx, nz=nz, triangles=tris,
            vertices=len(mesh.vertices), grid_cells=0, seconds=None,
            status=NOT_RUN, copies=copies, note=reason,
        )
    _prune_import_files(path)
    t0 = time.perf_counter()
    write_obj(mesh, path)
    elapsed = time.perf_counter() - t0
    size_mb = path.stat().st_size / (1024 * 1024)
    return SubTaskResult(
        subtask="obj_write", nx=nx, nz=nz, triangles=tris,
        vertices=len(mesh.vertices), grid_cells=0, copies=copies,
        seconds=round(elapsed, 6),
        throughput=round(tris / elapsed, 1) if elapsed > 0 else None,
        note=(
            f"file={size_mb:,.1f} MiB ({size_mb / elapsed:,.1f} MiB/s), "
            f"{size_mb * 1024 * 1024 / tris:.1f} B/tri, bodies={copies}"
        ),
        file_bytes=path.stat().st_size,
        file_mb_per_s=round(size_mb / elapsed, 3) if elapsed > 0 else None,
        peak_rss_mb=peak_rss_mb(),
        rss_mb=live_rss_mb(),
        mesh_mb=mesh_footprint_mb(mesh),
    )


def _run_obj_read(nx: int, nz: int, copies: int = 1) -> SubTaskResult:
    """Measure parsing an OBJ file back into a mesh, and check what came back.

    If the file is not on disk the sub-task prepares it first and records that
    preparation time in the note (untimed), so it can be run on its own.  The
    parsed mesh is checked against the mesh it came from -- vertex/face counts
    and the signed volume against the analytic reference -- because a fast parse
    that silently drops or duplicates triangles is a defect, not a measurement.
    """
    from tensorcad.export import read_obj, write_obj
    from tensorcad.hull import analytic_volume
    from tensorcad.mesh import total_signed_volume

    mesh = _mesh_for(nx, nz, copies)
    tris = len(mesh.faces)
    path = import_obj_path(nx, nz, copies)
    prepared = 0.0
    if not path.exists():
        reason = _obj_space_guard(tris)
        if reason is not None:
            return SubTaskResult(
                subtask="obj_read", nx=nx, nz=nz, triangles=tris,
                vertices=len(mesh.vertices), grid_cells=0, seconds=None,
                status=NOT_RUN, copies=copies,
                note=f"needs the OBJ file first; {reason}",
            )
        _prune_import_files(path)
        t_prep = time.perf_counter()
        write_obj(mesh, path)
        prepared = time.perf_counter() - t_prep

    t0 = time.perf_counter()
    parsed = read_obj(path)
    elapsed = time.perf_counter() - t0

    if len(parsed.vertices) != len(mesh.vertices) or len(parsed.faces) != tris:
        raise ValueError(
            f"OBJ round trip lost geometry: parsed {len(parsed.vertices)} vertices/"
            f"{len(parsed.faces)} faces, expected {len(mesh.vertices)}/{tris}"
        )

    volume = total_signed_volume(parsed)
    reference = copies * analytic_volume(10.0, 2.0, 0.8)
    rel_err = (volume - reference) / reference if reference else float("nan")
    size_mb = path.stat().st_size / (1024 * 1024)
    note = (
        f"file={size_mb:,.1f} MiB ({size_mb / elapsed:,.1f} MiB/s), "
        f"volume={volume:.6f} m^3, analytic={reference:.6f} m^3, "
        f"rel_err={rel_err:+.3e}, bodies={copies}"
    )
    if prepared:
        note += f", file prepared in {prepared:.3f}s (not counted)"
    return SubTaskResult(
        subtask="obj_read", nx=nx, nz=nz, triangles=len(parsed.faces),
        vertices=len(parsed.vertices), grid_cells=0, copies=copies,
        seconds=round(elapsed, 6),
        throughput=round(len(parsed.faces) / elapsed, 1) if elapsed > 0 else None,
        note=note,
        file_bytes=path.stat().st_size,
        file_mb_per_s=round(size_mb / elapsed, 3) if elapsed > 0 else None,
        volume_rel_err=rel_err,
        peak_rss_mb=peak_rss_mb(),
        rss_mb=live_rss_mb(),
    )


# ---------------------------------------------------------------------------
# Fibonacci directions
# ---------------------------------------------------------------------------

def _fibonacci_directions(n: int) -> list[tuple[float, float, float]]:
    """Deterministic unit vectors on a sphere (golden spiral)."""
    if n < 1:
        raise ValueError("need at least 1 direction")
    golden = (1.0 + math.sqrt(5.0)) / 2.0
    directions = []
    for i in range(n):
        theta = 2.0 * math.pi * i / golden
        phi = math.acos(1.0 - 2.0 * (i + 0.5) / n)
        st, ct = math.sin(theta), math.cos(theta)
        sp, cp = math.sin(phi), math.cos(phi)
        directions.append((ct * sp, st * sp, cp))
    return directions


# ---------------------------------------------------------------------------
# Benchmark runner
# ---------------------------------------------------------------------------

_SUBTASKS_MESH = {
    "load": _run_load,
    "point_locate": _run_point_locate,
    "point_locate_accel": _run_point_locate_accel,
    "volume": _run_volume,
    "obj_write": _run_obj_write,
    "obj_read": _run_obj_read,
}

#: Which sub-tasks accept the ray-query count argument.
_QUERY_SUBTASKS = ("point_locate", "point_locate_accel")

#: Above this triangle count the brute-force reference is skipped instead of
#: allowed to dominate the run (it is O(N) per ray: at 1M triangles a 10k-ray
#: sweep would take hours).
DEFAULT_MAX_NAIVE_TRIANGLES = 300_000


def run_perf_benchmark(
    ladder: tuple[tuple[int, int], ...] = DEFAULT_LADDER,
    voxel_ladder: tuple[int, ...] = VOXEL_LADDER,
    subtasks: tuple[str, ...] = (
        "load", "point_locate", "point_locate_accel", "voxelize", "volume",
        "obj_write", "obj_read",  # import/export sub-tasks (GAP-CAD-D4-PERF-SCALE §11)
    ),
    max_triangles: int | None = None,
    n_queries: int = QUICK_N_QUERIES,
    max_naive_triangles: int | None = DEFAULT_MAX_NAIVE_TRIANGLES,
) -> PerfBenchmarkResult:
    """Run the full performance benchmark ladder.

    Parameters
    ----------
    ladder:
        (nx, nz) pairs defining the mesh resolution ladder.
    voxel_ladder:
        Grid resolutions for the voxelization sub-task.
    subtasks:
        Which sub-tasks to run (subset of load/point_locate/voxelize/volume).
    max_triangles:
        Skip any step that would produce more than this many triangles.
    n_queries:
        Number of ray queries for point_locate / point_locate_accel.
    max_naive_triangles:
        Skip the brute-force ``point_locate`` reference above this triangle
        count (its cost grows as N x queries); ``point_locate_accel`` still
        covers such rungs.
    """
    t_start = time.perf_counter()
    results: list[SubTaskResult] = []
    cluster_warnings: list[str] = []

    # Mesh-based sub-tasks
    mesh_tasks = [t for t in subtasks if t != "voxelize"]
    for nx, nz in ladder:
        est_tris = 4 * (nx - 1) * (nz - 1)
        if max_triangles is not None and est_tris > max_triangles:
            continue

        for task_name in mesh_tasks:
            fn = _SUBTASKS_MESH.get(task_name)
            if fn is None:
                continue
            if (
                task_name == "point_locate"
                and max_naive_triangles is not None
                and est_tris > max_naive_triangles
            ):
                results.append(SubTaskResult(
                    subtask=task_name, nx=nx, nz=nz,
                    triangles=est_tris, vertices=0, grid_cells=0,
                    seconds=None, status=NOT_RUN,
                    note=(
                        f"brute-force reference skipped above "
                        f"{max_naive_triangles:,} triangles (O(N) per ray); "
                        f"use point_locate_accel"
                    ),
                ))
                continue
            try:
                if task_name in _QUERY_SUBTASKS:
                    result = fn(nx, nz, n_queries=n_queries)
                else:
                    result = fn(nx, nz)
                results.append(result)
                if result.triangles >= CLUSTER_THRESHOLD:
                    cluster_warnings.append(
                        f"{task_name} at {nx}x{nz} ({result.triangles:,} triangles) "
                        f"took {result.seconds:.3f}s — consider cluster for larger scales"
                    )
            except Exception as exc:
                results.append(SubTaskResult(
                    subtask=task_name, nx=nx, nz=nz,
                    triangles=0, vertices=0, grid_cells=0,
                    seconds=None, status=FAILED, note=str(exc),
                ))

    # Voxelization sub-task
    if "voxelize" in subtasks:
        for grid_res in voxel_ladder:
            cells = grid_res ** 3
            if max_triangles is not None and cells > max_triangles:
                continue
            try:
                result = _run_voxelize(grid_res)
                results.append(result)
                if cells >= CLUSTER_THRESHOLD:
                    cluster_warnings.append(
                        f"voxelize at {grid_res}^3 ({cells:,} cells) "
                        f"took {result.seconds:.3f}s — consider cluster for larger scales"
                    )
            except Exception as exc:
                results.append(SubTaskResult(
                    subtask="voxelize", nx=grid_res, nz=0,
                    triangles=0, vertices=0, grid_cells=0,
                    seconds=None, status=FAILED, note=str(exc),
                ))

    elapsed = time.perf_counter() - t_start
    return PerfBenchmarkResult(
        results=tuple(results),
        gates=scaling_gates(results),
        environment=_environment(),
        seconds=round(elapsed, 3),
        cluster_warnings=tuple(cluster_warnings),
    )


# ---------------------------------------------------------------------------
# 10^7-band fleet runner
# ---------------------------------------------------------------------------

#: Default base rung for the fleet tier: one 4.0e6-triangle generated hull,
#: assembled into 2 (8.0e6) and 3 (1.2e7) bodies.  Two fleet sizes from one
#: generation give the 10^6 -> 10^7 doubling without paying for a second
#: multi-minute generation.
FLEET_BASE_LADDER: tuple[tuple[int, int], ...] = ((2001, 501),)
FLEET_COPIES_LADDER: tuple[int, ...] = (2, 3)

#: Sub-tasks that are meaningful on a multi-body fleet rung.  ``point_locate``
#: (brute-force star-shaped reference) and ``voxelize`` (a grid-resolution
#: ladder, not a mesh one) are structurally inapplicable: running them would
#: put rows in the table and JSON whose status is ``failed``, which reads as a
#: measurement error rather than "does not apply".  They are recorded once per
#: base rung with ``status=not_run`` and a reason instead.
FLEET_SUBTASKS: tuple[str, ...] = (
    "point_locate_accel", "volume", "obj_write", "obj_read",
)


def run_fleet_benchmark(
    base_ladder: tuple[tuple[int, int], ...] = FLEET_BASE_LADDER,
    copies_ladder: tuple[int, ...] = FLEET_COPIES_LADDER,
    subtasks: tuple[str, ...] = ("load", "fleet_assemble", "point_locate_accel", "volume"),
    n_queries: int = 2000,
    max_triangles: int | None = None,
) -> PerfBenchmarkResult:
    """Measure the 10^6-10^7 band on multi-body (fleet) meshes.

    One hull is generated per base rung (``load``, timed); each body count then
    gets ``fleet_assemble`` (the O(N) copy), ``point_locate_accel`` and
    ``volume``.  The brute-force ``point_locate`` reference is deliberately not
    offered here: ``_RayCaster`` assumes a star-shaped single body.
    """
    t_start = time.perf_counter()
    results: list[SubTaskResult] = []
    cluster_warnings: list[str] = []
    reported_skips: set[tuple[str, str]] = set()

    for nx, nz in base_ladder:
        base_tris = 4 * (nx - 1) * (nz - 1)
        if max_triangles is not None and base_tris * max(copies_ladder) > max_triangles:
            continue
        try:
            load = _run_load(nx, nz)
        except Exception as exc:  # pragma: no cover - defensive
            results.append(SubTaskResult(
                subtask="load", nx=nx, nz=nz, triangles=0, vertices=0,
                grid_cells=0, seconds=None, status=FAILED, note=str(exc),
            ))
            continue
        results.append(load)
        for copies in copies_ladder:
            try:
                results.append(_run_fleet_assemble(nx, nz, copies))
            except Exception as exc:
                results.append(SubTaskResult(
                    subtask="fleet_assemble", nx=nx, nz=nz, triangles=0, vertices=0,
                    grid_cells=0, seconds=None, status=FAILED, note=str(exc), copies=copies,
                ))
                continue
            for task_name in subtasks:
                if task_name in ("load", "fleet_assemble"):
                    continue
                if task_name not in FLEET_SUBTASKS:
                    key = ("skip", task_name)
                    if key not in reported_skips:
                        reported_skips.add(key)
                        results.append(SubTaskResult(
                            subtask=task_name, nx=nx, nz=nz, triangles=0, vertices=0,
                            grid_cells=0, seconds=None, status=NOT_RUN, copies=copies,
                            note=(
                                "not applicable to a multi-body fleet rung: "
                                + ("the brute-force reference assumes a star-shaped single body"
                                   if task_name == "point_locate" else
                                   "the sub-task is a grid-resolution ladder, not a mesh one")
                            ),
                        ))
                    continue
                try:
                    if task_name == "point_locate_accel":
                        r = _run_point_locate_accel(nx, nz, n_queries=n_queries, copies=copies)
                    elif task_name == "volume":
                        r = _run_volume(nx, nz, copies=copies)
                    else:
                        r = _SUBTASKS_MESH[task_name](nx, nz, copies=copies)
                    results.append(r)
                    if r.triangles >= CLUSTER_THRESHOLD:
                        cluster_warnings.append(
                            f"{task_name} at {nx}x{nz} x{copies} bodies "
                            f"({r.triangles:,} triangles) took {r.seconds:.3f}s, "
                            f"peak RSS {r.peak_rss_mb or 0:.0f} MiB"
                        )
                except Exception as exc:
                    results.append(SubTaskResult(
                        subtask=task_name, nx=nx, nz=nz, triangles=0, vertices=0,
                        grid_cells=0, seconds=None, status=FAILED, note=str(exc),
                        copies=copies,
                    ))
        clear_mesh_cache()
        # The fleet cache is dropped per rung too: one 10^7 mesh is several GB.
        clear_fleet_cache()

    elapsed = time.perf_counter() - t_start
    return PerfBenchmarkResult(
        results=tuple(results),
        gates=scaling_gates(results),
        environment=_environment(),
        seconds=round(elapsed, 3),
        cluster_warnings=tuple(cluster_warnings),
    )


# ---------------------------------------------------------------------------
# Regression gates
# ---------------------------------------------------------------------------

#: Prefix of the scale-derived (complexity) gates, so they can be recognised and
#: preserved when wall-clock regression gates are added on top.
SCALING_GATE_PREFIX = "scaling_"

#: Upper bound on the fitted log-log exponent, per sub-task.  A wall-clock gate
#: against a baseline measured on tiny meshes cannot see an algorithmic
#: regression; the exponent can.  Measured exponents on the committed ladder are
#: ~1.0 for the O(N) paths (``load``, ``volume``, brute-force ``point_locate``)
#: and ~0.3 for the indexed point location; these limits leave head-room for
#: shared-machine noise while still flagging a regression towards O(N^2)
#: (exponent ~2).
SCALING_LIMITS: dict[str, float] = {
    "load": 1.5,
    "fleet_assemble": 1.5,
    "point_locate": 1.5,
    "point_locate_accel": 0.8,
    "volume": 1.5,
    "voxelize": 1.5,
    "obj_write": 1.5,
    "obj_read": 1.5,
}

#: Fewest distinct sizes needed before an exponent is gated at all.
MIN_SCALING_POINTS = 3


def measurement_key(subtask: str, nx: int, nz: int, copies: int = 1) -> str:
    """Stable identity of one measurement, used as the baseline/gate key.

    Single-body keys keep the historical ``subtask:nxXxnz`` spelling so existing
    baselines stay valid; fleet rungs append the body count.
    """
    key = f"{subtask}:{nx}x{nz}"
    if copies != 1:
        key += f":b{copies}"
    return key


def _scale_axis(result: SubTaskResult) -> float:
    """The size axis of a measurement: triangles for meshes, cells for voxels."""
    if result.subtask == "voxelize":
        return float(result.grid_cells)
    return float(result.triangles)


def _loglog_slope(points: list[tuple[float, float]]) -> tuple[float, float]:
    """Least-squares slope of log(seconds) vs log(size), and the R^2 of the fit."""
    xs = [math.log(x) for x, _ in points]
    ys = [math.log(y) for _, y in points]
    n = len(points)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0.0:
        raise ValueError("degenerate size axis (all sizes equal)")
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = my - slope * mx
    ss_res = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - my) ** 2 for y in ys)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0
    return slope, r_squared


def scaling_gates(
    results: tuple[SubTaskResult, ...] | list[SubTaskResult],
    limits: dict[str, float] | None = None,
    min_points: int = MIN_SCALING_POINTS,
) -> tuple[GateVerdict, ...]:
    """Gate the measured complexity exponent of every sub-task.

    This is the *scale-invariant* half of the performance gate.  The wall-clock
    gate (:func:`evaluate_gates`) compares one run against a stored baseline at
    the same size, which is why its committed baseline only covers the cheap
    ladder: a 10^6-triangle run is too slow for CI, so that band had no gate at
    all.  Fitting the exponent closes that hole cheaply - a change that turns an
    O(N) pass into O(N^2) shows up as exponent ~2 even on small rungs.

    A sub-task with fewer than ``min_points`` distinct sizes is reported with
    ``status=NOT_RUN`` (explicitly "not checked", never a silent pass).
    """
    table = SCALING_LIMITS if limits is None else limits
    groups: dict[tuple[str, int], dict[float, float]] = {}
    for r in results:
        if r.seconds is None or r.seconds <= 0.0:
            continue
        size = _scale_axis(r)
        if size <= 0.0:
            continue
        # Keep the slowest observation per exact size (repeat runs).
        bucket = groups.setdefault((r.subtask, r.copies), {})
        bucket[size] = max(bucket.get(size, 0.0), r.seconds)

    gates: list[GateVerdict] = []
    for (subtask, copies), bucket in sorted(groups.items()):
        suffix = "" if copies == 1 else f"[x{copies}]"
        name = f"{SCALING_GATE_PREFIX}{subtask}{suffix}"
        limit = table.get(subtask)
        points = sorted(bucket.items())
        if limit is None:
            gates.append(GateVerdict(
                name=name, passed=True, status=NOT_RUN,
                detail=f"not_run: no exponent limit configured for {subtask!r}",
            ))
            continue
        if len(points) < min_points:
            gates.append(GateVerdict(
                name=name, passed=True, status=NOT_RUN,
                detail=(
                    f"not_run: {len(points)} distinct size(s), need {min_points} "
                    f"to fit an exponent"
                ),
            ))
            continue
        slope, r_squared = _loglog_slope(points)
        passed = slope <= limit
        gates.append(GateVerdict(
            name=name,
            passed=passed,
            detail=(
                f"exponent={slope:.2f} (limit {limit:.2f}), R^2={r_squared:.3f}, "
                f"{len(points)} sizes {int(points[0][0]):,}..{int(points[-1][0]):,}"
            ),
            ratio=round(slope, 3),
        ))
    return tuple(gates)


#: Environment fields that make two wall-clock measurements incomparable.  A
#: stored baseline is a *time*, and time only means something on the machine it
#: was taken on, so a baseline recorded on another architecture must not be
#: allowed to flag RED (or, worse, to pass) on this one.
ENV_HARD_KEYS = ("machine", "platform")
#: Environment fields that bias a comparison but do not invalidate it; they are
#: surfaced in the gate detail instead of silently skewing the ratio.
ENV_SOFT_KEYS = ("cpu_count", "python")


def environment_compat(
    baseline_env: dict[str, Any] | None,
    current_env: dict[str, Any] | None,
) -> tuple[bool, str]:
    """Can wall-clock times from ``baseline_env`` be compared with ``current_env``?

    Returns ``(comparable, reason)``.  A baseline without an environment record
    (older baselines) is assumed comparable and said so in ``reason``; a
    different architecture / OS is a hard mismatch; core count and interpreter
    version only annotate the result.
    """
    if not baseline_env or not current_env:
        return True, "no environment record on one side (assumed comparable)"
    for key in ENV_HARD_KEYS:
        b, c = baseline_env.get(key), current_env.get(key)
        if b and c and b != c:
            return False, f"{key} differs: measured {c!r} vs baseline {b!r}"
    notes = []
    for key in ENV_SOFT_KEYS:
        b, c = baseline_env.get(key), current_env.get(key)
        if b and c and b != c:
            notes.append(f"{key} {c} vs baseline {b}")
    return True, "; ".join(notes) if notes else "environment matches baseline"


def evaluate_gates(
    result: PerfBenchmarkResult,
    baseline_path: str | Path | None = None,
    regression_factor: float = 2.0,
) -> tuple[PerfBenchmarkResult, list[str]]:
    """Add regression gates to a benchmark result.

    If ``baseline_path`` is given, each sub-task at each (nx, nz) is compared
    against the baseline time.  Exceeding ``regression_factor`` × baseline
    produces a RED gate.
    """
    baseline: dict[str, float] = {}
    baseline_env: dict[str, Any] | None = None
    if baseline_path is not None:
        bp = Path(baseline_path)
        if bp.exists():
            data = json.loads(bp.read_text(encoding="utf-8"))
            baseline_env = data.get("environment")
            for entry in data.get("results", []):
                key = measurement_key(
                    entry["subtask"], entry["nx"], entry["nz"], entry.get("copies", 1)
                )
                if entry.get("seconds") is not None:
                    baseline[key] = entry["seconds"]

    gates: list[GateVerdict] = []
    comparable, env_reason = environment_compat(baseline_env, result.environment)
    if baseline and not comparable:
        # Refuse the comparison instead of inventing a RED/pass number: a
        # wall-clock baseline is machine-specific.  status=not_run (never a
        # pass) so the verdict cannot be reported as green.
        scaling_only = tuple(g for g in result.gates if g.name.startswith(SCALING_GATE_PREFIX))
        refusal = GateVerdict(
            name="regression_environment",
            passed=True,
            status=NOT_RUN,
            detail=f"not_run: baseline not comparable ({env_reason})",
        )
        return replace(result, gates=scaling_only + (refusal,)), []
    for r in result.results:
        if r.seconds is None:
            continue
        key = measurement_key(r.subtask, r.nx, r.nz, r.copies)
        if key not in baseline:
            continue
        bl = baseline[key]
        ratio = r.seconds / bl if bl > 0 else float("inf")
        passed = ratio <= regression_factor
        gates.append(GateVerdict(
            name=f"regression_{key}",
            passed=passed,
            detail=(
                f"{r.seconds:.4f}s vs {bl:.4f}s baseline (ratio={ratio:.2f}x, "
                f"limit={regression_factor}x; env: {env_reason})"
            ),
            baseline_seconds=bl,
            measured_seconds=r.seconds,
            ratio=round(ratio, 3),
        ))

    # Preserve the scale-derived gates computed by the runner; this function
    # only appends the wall-clock comparison against a stored baseline.
    scaling = tuple(g for g in result.gates if g.name.startswith(SCALING_GATE_PREFIX))
    gates = list(scaling) + gates
    return replace(result, gates=tuple(gates)), [
        g.name for g in gates if not g.passed
    ]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="TensorCAD large-scale geometry performance benchmark")
    parser.add_argument("--quick", action="store_true", help="Use quick ladder (up to ~120k triangles)")
    parser.add_argument("--full", action="store_true", help="Use full ladder (up to ~18M triangles)")
    parser.add_argument("--ladder", type=str, help="Custom ladder as 'nx1,nz1;nx2,nz2;...'")
    parser.add_argument("--subtasks", type=str,
                        default="load,point_locate,point_locate_accel,voxelize,volume",
                        help="Comma-separated subtask names "
                             "(load,point_locate,point_locate_accel,voxelize,volume,obj_write,obj_read)")
    parser.add_argument("--max-triangles", type=int, default=None,
                        help="Skip steps producing more than this many triangles")
    parser.add_argument("--n-queries", type=int, default=None,
                        help="Number of ray queries for point_locate")
    parser.add_argument("--fleet-tier", action="store_true",
                        help="Measure the 10^6-10^7 band on multi-body fleet meshes "
                             "(see build_fleet_mesh); uses --ladder/--fleet-copies")
    parser.add_argument("--fleet-copies", type=str, default="2,3",
                        help="Comma-separated body counts for --fleet-tier (default 2,3)")
    parser.add_argument("--max-naive-triangles", type=int, default=DEFAULT_MAX_NAIVE_TRIANGLES,
                        help="Skip the brute-force point_locate reference above this triangle count")
    parser.add_argument("--save-baseline", type=str, default=None,
                        help="Save results as baseline JSON at this path")
    parser.add_argument("--baseline", type=str, default=None,
                        help="Compare against baseline JSON at this path")
    parser.add_argument("--regression-factor", type=float, default=DEFAULT_GATE_FACTOR,
                        help="RED threshold: measured > factor × baseline")
    parser.add_argument("--output", type=str, default=None,
                        help="Write result JSON to this path")
    args = parser.parse_args(argv)

    if args.quick:
        ladder = QUICK_LADDER
        voxel_ladder = QUICK_VOXEL_LADDER
        n_queries = args.n_queries or QUICK_N_QUERIES
    elif args.full:
        ladder = DEFAULT_LADDER
        voxel_ladder = VOXEL_LADDER
        n_queries = args.n_queries or FULL_N_QUERIES
    elif args.ladder:
        ladder = tuple(
            (int(a), int(b))
            for pair in args.ladder.split(";")
            for a, b in [pair.split(",")]
        )
        voxel_ladder = VOXEL_LADDER
        n_queries = args.n_queries or QUICK_N_QUERIES
    else:
        ladder = QUICK_LADDER
        voxel_ladder = QUICK_VOXEL_LADDER
        n_queries = QUICK_N_QUERIES

    subtasks = tuple(s.strip() for s in args.subtasks.split(","))

    if args.fleet_tier:
        copies_ladder = tuple(int(c) for c in args.fleet_copies.split(",") if c.strip())
        if args.ladder and not (args.quick or args.full):
            base_ladder = ladder
        else:
            base_ladder = FLEET_BASE_LADDER
        fleet_subtasks = tuple(t for t in subtasks if t in ("load",) + FLEET_SUBTASKS) or ("load",) + FLEET_SUBTASKS
        result = run_fleet_benchmark(
            base_ladder=base_ladder,
            copies_ladder=copies_ladder,
            subtasks=fleet_subtasks,
            n_queries=n_queries or FULL_N_QUERIES,
            max_triangles=args.max_triangles,
        )
    else:
        result = run_perf_benchmark(
            ladder=ladder,
            voxel_ladder=voxel_ladder,
            subtasks=subtasks,
            max_triangles=args.max_triangles,
            n_queries=n_queries,
            max_naive_triangles=args.max_naive_triangles,
        )

    failed_gates = []
    if args.baseline:
        result, failed_gates = evaluate_gates(result, args.baseline, args.regression_factor)

    if args.save_baseline:
        result.write_json(args.save_baseline)
        print(f"Baseline saved to {args.save_baseline}")

    if args.output:
        result.write_json(args.output)
        print(f"Results written to {args.output}")

    # Print summary table
    _print_summary(result)

    if failed_gates:
        print(f"\nRED gates ({len(failed_gates)}):")
        for name in failed_gates:
            print(f"  ✗ {name}")
        return 1

    if not result.gates:
        print("\nNo baseline for regression gates.")
    elif result.verdict == "pass":
        extra = f" ({result.gates_not_run} not_run, not counted as a pass)" if result.gates_not_run else ""
        print(f"\nAll gates passed.{extra}")
    else:
        print(
            f"\nNo gate could be checked: all {result.gates_not_run} gate(s) not_run. "
            f"This is NOT a pass."
        )
    return 0


def _print_summary(result: PerfBenchmarkResult) -> None:
    """Print a human-readable summary table."""
    print(f"\n{'='*80}")
    print(f"TensorCAD Performance Benchmark  (schema {BENCHMARK_SCHEMA})")
    print(f"{'='*80}")
    print(f"Total time: {result.seconds:.3f}s")
    print(f"Environment: Python {result.environment.get('python', '?')}, "
          f"TensorCAD {result.environment.get('tensorcad', '?')}, "
          f"CPU={result.environment.get('cpu_count', '?')}")
    la = result.environment.get("loadavg", {})
    if la:
        print(f"Load avg: {la.get('load1', '?')}/{la.get('load5', '?')}/{la.get('load15', '?')}")
    print()

    # Group by subtask
    by_task: dict[str, list[SubTaskResult]] = {}
    for r in result.results:
        by_task.setdefault(r.subtask, []).append(r)

    for task_name in ("load", "fleet_assemble", "point_locate", "point_locate_accel",
                      "voxelize", "volume", "obj_write", "obj_read"):
        entries = by_task.get(task_name, [])
        if not entries:
            continue
        print(f"--- {task_name} ---")
        if task_name == "voxelize":
            print(f"  {'grid':>6} {'cells':>12} {'seconds':>10} {'throughput':>14} {'status':>10}")
            for r in entries:
                tp = f"{r.throughput:,.0f}" if r.throughput else "N/A"
                sec = f"{r.seconds:.4f}" if r.seconds is not None else "N/A"
                print(f"  {r.nx:>6} {r.grid_cells:>12,} {sec:>10} {tp:>14} {r.status:>10}")
        else:
            print(f"  {'nx':>5} {'nz':>5} {'bodies':>6} {'triangles':>11} {'vertices':>11} {'seconds':>10} {'throughput':>14} {'build':>8} {'rss_MB':>8} {'mesh_MB':>8} {'status':>10}")
            for r in entries:
                tp = f"{r.throughput:,.0f}" if r.throughput else "N/A"
                sec = f"{r.seconds:.4f}" if r.seconds is not None else "N/A"
                bld = f"{r.build_seconds:.3f}" if r.build_seconds is not None else "-"
                rss = f"{r.peak_rss_mb:,.0f}" if r.peak_rss_mb is not None else "-"
                meshm = f"{r.mesh_mb:,.0f}" if r.mesh_mb is not None else "-"
                print(f"  {r.nx:>5} {r.nz:>5} {r.copies:>6} {r.triangles:>11,} {r.vertices:>11,} {sec:>10} {tp:>14} {bld:>8} {rss:>8} {meshm:>8} {r.status:>10}")
        print()

    if result.cluster_warnings:
        print("⚠ Cluster warnings:")
        for w in result.cluster_warnings:
            print(f"  {w}")
        print()

    if result.gates:
        print("--- Regression Gates (wall-clock ratio + measured complexity exponent) ---")
        for g in result.gates:
            if g.status == NOT_RUN:
                mark = "· not_run"
            else:
                mark = "✓" if g.passed else "✗ RED"
            print(f"  {mark:>9} {g.name}: {g.detail}")
        print()


if __name__ == "__main__":
    sys.exit(main())
