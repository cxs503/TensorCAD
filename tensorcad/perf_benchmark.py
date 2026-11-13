"""Large-scale geometry performance benchmark (10^5--10^7 triangle regime).

Answers the question in GAP-CAD-D4-PERF-SCALE: *at what mesh size does each
geometry operation become expensive, and when does it need a cluster?*

Four sub-tasks are measured on the same mesh ladder:

``load``          Mesh construction from analytic parameters (hull generation).
``point_locate``  Ray-casting / signed-distance queries (point-in-geometry).
``voxelize``      SDF sampling onto a regular grid (to_voxel_grid).
``volume``        Surface-integral volume by the divergence theorem.

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

    # Save baseline for future regression gates
    python -m tensorcad.perf_benchmark --full --save-baseline runs/perf_baseline.json

    # Check against a saved baseline
    python -m tensorcad.perf_benchmark --full --baseline runs/perf_baseline.json

Exit code is 0 when all gates pass, 1 when any gate fails.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Data contracts
# ---------------------------------------------------------------------------

BENCHMARK_SCHEMA = "tensorcad.perf_benchmark/1"

#: Status tags for individual measurements.
MEASURED = "measured"
NOT_RUN = "not_run"
FAILED = "failed"

#: The mesh ladder: (nx, nz) pairs that produce increasing triangle counts.
#: The default ladder goes from ~2k to ~10M triangles.
DEFAULT_LADDER: tuple[tuple[int, int], ...] = (
    (21, 9),      # ~2k triangles
    (41, 13),     # ~8k triangles
    (81, 21),     # ~30k triangles
    (161, 41),    # ~120k triangles
    (321, 81),    # ~480k triangles
    (641, 161),   # ~1.9M triangles
    (1281, 321),  # ~7.6M triangles
    (2001, 501),  # ~18M triangles (upper bound, may be skipped)
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
CLUSTER_THRESHOLD = 5_000_000

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

    def to_dict(self) -> dict[str, Any]:
        return {
            "subtask": self.subtask,
            "nx": self.nx,
            "nz": self.nz,
            "triangles": self.triangles,
            "vertices": self.vertices,
            "grid_cells": self.grid_cells,
            "seconds": self.seconds,
            "status": self.status,
            "throughput": self.throughput,
            "note": self.note,
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

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "passed": self.passed,
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

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": BENCHMARK_SCHEMA,
            "passed": self.passed,
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
# Sub-task implementations
# ---------------------------------------------------------------------------

def _run_load(nx: int, nz: int) -> SubTaskResult:
    """Measure hull generation time."""
    from tensorcad.hull import HullParameters, generate_hull
    params = HullParameters(length=10.0, beam=2.0, draft=0.8, nx=nx, nz=nz)
    t0 = time.perf_counter()
    model = generate_hull(params)
    elapsed = time.perf_counter() - t0
    tris = len(model.mesh.faces)
    verts = len(model.mesh.vertices)
    return SubTaskResult(
        subtask="load", nx=nx, nz=nz, triangles=tris, vertices=verts,
        grid_cells=0,
        seconds=round(elapsed, 6),
        throughput=round(tris / elapsed, 1) if elapsed > 0 else None,
    )


def _run_point_locate(nx: int, nz: int, n_queries: int = 1000) -> SubTaskResult:
    """Measure ray-casting / point-location throughput.

    The _RayCaster casts rays from the origin; the hull is centered at the
    origin so rays in all directions that intersect the hull will return a
    radius.  We use a Fibonacci sphere for uniform directional coverage.
    """
    from tensorcad.hull import HullParameters, generate_hull
    from tensorcad.benchmark import _RayCaster

    params = HullParameters(length=10.0, beam=2.0, draft=0.8, nx=nx, nz=nz)
    model = generate_hull(params)
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
        note=f"hits={hits}/{n_queries}",
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
    )


def _run_volume(nx: int, nz: int) -> SubTaskResult:
    """Measure volume integration time."""
    from tensorcad.hull import HullParameters, generate_hull
    from tensorcad import mesh as meshlib

    params = HullParameters(length=10.0, beam=2.0, draft=0.8, nx=nx, nz=nz)
    model = generate_hull(params)
    tris = len(model.mesh.faces)
    verts = len(model.mesh.vertices)

    t0 = time.perf_counter()
    vol = meshlib.total_signed_volume(model.mesh)
    elapsed = time.perf_counter() - t0

    return SubTaskResult(
        subtask="volume", nx=nx, nz=nz, triangles=tris, vertices=verts,
        grid_cells=0,
        seconds=round(elapsed, 6),
        throughput=round(tris / elapsed, 1) if elapsed > 0 else None,
        note=f"volume={vol:.6f} m^3",
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
    "volume": _run_volume,
}


def run_perf_benchmark(
    ladder: tuple[tuple[int, int], ...] = DEFAULT_LADDER,
    voxel_ladder: tuple[int, ...] = VOXEL_LADDER,
    subtasks: tuple[str, ...] = ("load", "point_locate", "voxelize", "volume"),
    max_triangles: int | None = None,
    n_queries: int = QUICK_N_QUERIES,
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
        Number of ray queries for point_locate.
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
            try:
                if task_name == "point_locate":
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
        environment=_environment(),
        seconds=round(elapsed, 3),
        cluster_warnings=tuple(cluster_warnings),
    )


# ---------------------------------------------------------------------------
# Regression gate
# ---------------------------------------------------------------------------

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
    if baseline_path is not None:
        bp = Path(baseline_path)
        if bp.exists():
            data = json.loads(bp.read_text(encoding="utf-8"))
            for entry in data.get("results", []):
                key = f"{entry['subtask']}:{entry['nx']}x{entry['nz']}"
                if entry.get("seconds") is not None:
                    baseline[key] = entry["seconds"]

    gates: list[GateVerdict] = []
    for r in result.results:
        if r.seconds is None:
            continue
        key = f"{r.subtask}:{r.nx}x{r.nz}"
        if key not in baseline:
            continue
        bl = baseline[key]
        ratio = r.seconds / bl if bl > 0 else float("inf")
        passed = ratio <= regression_factor
        gates.append(GateVerdict(
            name=f"regression_{key}",
            passed=passed,
            detail=f"{r.seconds:.4f}s vs {bl:.4f}s baseline (ratio={ratio:.2f}x, limit={regression_factor}x)",
            baseline_seconds=bl,
            measured_seconds=r.seconds,
            ratio=round(ratio, 3),
        ))

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
    parser.add_argument("--subtasks", type=str, default="load,point_locate,voxelize,volume",
                        help="Comma-separated subtask names")
    parser.add_argument("--max-triangles", type=int, default=None,
                        help="Skip steps producing more than this many triangles")
    parser.add_argument("--n-queries", type=int, default=None,
                        help="Number of ray queries for point_locate")
    parser.add_argument("--save-baseline", type=str, default=None,
                        help="Save results as baseline JSON at this path")
    parser.add_argument("--baseline", type=str, default=None,
                        help="Compare against baseline JSON at this path")
    parser.add_argument("--regression-factor", type=float, default=2.0,
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

    result = run_perf_benchmark(
        ladder=ladder,
        voxel_ladder=voxel_ladder,
        subtasks=subtasks,
        max_triangles=args.max_triangles,
        n_queries=n_queries,
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

    print("\nAll gates passed." if result.gates else "\nNo baseline for regression gates.")
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

    for task_name in ("load", "point_locate", "voxelize", "volume"):
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
            print(f"  {'nx':>5} {'nz':>5} {'triangles':>10} {'vertices':>10} {'seconds':>10} {'throughput':>14} {'status':>10}")
            for r in entries:
                tp = f"{r.throughput:,.0f}" if r.throughput else "N/A"
                sec = f"{r.seconds:.4f}" if r.seconds is not None else "N/A"
                print(f"  {r.nx:>5} {r.nz:>5} {r.triangles:>10,} {r.vertices:>10,} {sec:>10} {tp:>14} {r.status:>10}")
        print()

    if result.cluster_warnings:
        print("⚠ Cluster warnings:")
        for w in result.cluster_warnings:
            print(f"  {w}")
        print()

    if result.gates:
        print("--- Regression Gates ---")
        for g in result.gates:
            mark = "✓" if g.passed else "✗ RED"
            print(f"  {mark} {g.name}: {g.detail}")
        print()


if __name__ == "__main__":
    sys.exit(main())
