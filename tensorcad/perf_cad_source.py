"""Performance ladder measured on a **real CAD file**, not on a generator.

Why this module exists
----------------------
``tensorcad.perf_benchmark`` answers "what do 10^6--10^7 triangles cost?" for
geometry this repository *generates* (``generate_hull``, subdivision, fleets of
copies).  Every table in ``docs/perf_scale_*.md`` therefore carries the same
caveat, stated in ``docs/perf_scale_import_10m.md`` §5.2: the band has never
been reached from an actual IGES/STEP file, and *real trimmed-surface
tessellation with tolerance welding and chunked parsing* is the last structural
hole of anchor ``GAP-CAD-D4-PERF-SCALE``.  This module closes it by running the
same style of measurement on a CAD source the repository already tracks:

``work/assets/onrt.igs``
    ONR tumblehome hull, 132 B-spline surfaces / 80 trimmed patches / 447 sewn
    edges, 3.4 MiB.  Patch-set route (untrimmed 128 patches, sewn): the triangle
    count is driven by the surface sampling density ``(u_samples, v_samples)``,
    which is the axis a CAD importer exposes.

``benchmarks/anchors/sources/KP505_modelScale_250mm.igs``
    KP505 propeller, 4 surfaces / 8 trimmed patches.  Trimmed route: the
    triangle count is driven by the **chord tolerance** (metres), i.e. the axis a
    real tessellator exposes.  It is the fidelity axis: a tolerance means
    nothing until you know how many triangles it costs.

Sub-tasks (all measured on the tessellated CAD mesh)
----------------------------------------------------
``cad_parse``       read the file into a B-rep / patch set.
``cad_tessellate``  sample the surfaces into triangles.
``cad_tessellate_chunked``
                    the same tessellation with the per-face sample grids
                    streamed in chunks (:meth:`BRep.tessellate`'s
                    ``chunk_faces``); the mesh must come back **bit-identical**,
                    so the only thing chunking buys is a bounded sampling peak.
``cad_weld``        tolerance welding (``Mesh.weld``): how many vertices merge,
                    how many open boundary edges close, whether the result is
                    watertight.
``cad_volume``      divergence-theorem volume of the tessellated mesh.
``cad_locate``      accelerated point location (uniform grid + DDA rays).

Gates
-----
* **wall-clock regression, 2x = RED** - per ``(subtask, source digest, route,
  parameter)`` key against a stored baseline.  A rung with no baseline key is
  ``not_run``, never a pass.
* **scaling exponent** - delegated to
  :func:`tensorcad.perf_benchmark.scaling_gates`, so the complexity half of the
  gate is identical to the generated-geometry ladder.

Everything is stdlib + this repository.  Run it as::

    python -m tensorcad.perf_cad_source --quick
    python -m tensorcad.perf_cad_source --source work/assets/onrt.igs \
        --rungs 8x6,16x12,32x24,64x48 --save-baseline benchmarks/perf_scale/cad_source_baseline.json
    python -m tensorcad.perf_cad_source --source work/assets/onrt.igs \
        --rungs 8x6,16x12 --baseline benchmarks/perf_scale/cad_source_baseline.json

Exit code 0 when every evaluated gate passes, 1 when any gate is RED.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

from .perf_benchmark import (
    DEFAULT_GATE_FACTOR,
    MEASURED,
    NOT_RUN,
    GateVerdict,
    SubTaskResult,
    live_rss_mb,
    mesh_footprint_mb,
    peak_rss_mb,
    scaling_gates,
)

SCHEMA = "tensorcad.perf_cad_source/1"

#: The two axes a CAD importer exposes (and therefore the two ladders).
ROUTE_PATCH = "patch"
ROUTE_TRIMMED = "trimmed"
ROUTES = (ROUTE_PATCH, ROUTE_TRIMMED)

#: Default sources, both tracked in this repository.
DEFAULT_PATCH_SOURCE = "work/assets/onrt.igs"
DEFAULT_TRIMMED_SOURCE = "benchmarks/anchors/sources/KP505_modelScale_250mm.igs"

#: Chord tolerances (metres) for the trimmed ladder of the 250 mm propeller.
DEFAULT_CHORD_TOLERANCES = (2e-2, 1e-2, 5e-3, 1e-3)

#: Sampling densities for the patch-set ladder.
DEFAULT_SAMPLE_RUNGS = ((8, 6), (16, 12), (32, 24), (64, 48))

#: A rung whose tessellation is expected to exceed this is not attempted (the
#: rung is recorded as ``not_run`` with the expected triangle count, never
#: silently skipped and never run into an OOM kill).
DEFAULT_MAX_TRIANGLES = 16_000_000

#: Vertices welded by exact position count as already shared; the tolerance
#: weld below uses the same default the import pipeline uses (1e-6 * extent).
DEFAULT_WELD_RELATIVE_TOLERANCE = 1e-6

#: Ray count for the point-location sub-task, matching the generated ladder.
DEFAULT_QUERIES = 1000

#: Faces per chunk for the chunked tessellation sub-task.
DEFAULT_CHUNK_FACES = 8


# ---------------------------------------------------------------------------
# source identity
# ---------------------------------------------------------------------------
def source_digest(path: str | Path) -> str:
    """SHA-256 of the source file: the identity every gate key hangs off."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def measurement_key(subtask: str, digest: str, route: str, parameter: str) -> str:
    """Stable key for the wall-clock baseline.

    The parameter is part of the key on purpose: A 2x gate that merged
    ``u=16,v=12`` with ``u=64,v=48`` would compare two different workloads and
    fire on nothing but the workload change.
    """
    return f"{subtask}|{digest[:16]}|{route}|{parameter}"


# ---------------------------------------------------------------------------
# rungs
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Rung:
    """One step of a ladder.

    ``u_samples``/``v_samples`` drive the patch-set route; ``chord_tolerance``
    (metres) drives the trimmed route.  Which field is meaningful is decided by
    the route, and :meth:`parameter` is what the gate keys on.
    """

    u_samples: int = 0
    v_samples: int = 0
    chord_tolerance: float | None = None

    def parameter(self, route: str) -> str:
        if route == ROUTE_TRIMMED:
            if self.chord_tolerance is None:
                raise ValueError("trimmed rung needs chord_tolerance")
            return f"chord={self.chord_tolerance:g}"
        return f"uv={self.u_samples}x{self.v_samples}"

    def expected_triangles(self, face_count: int) -> int:
        """Upper bound on the triangle count (2 quads per cell, trimmed <= whole)."""
        return 2 * self.u_samples * self.v_samples * face_count


def parse_rungs(spec: str) -> tuple[Rung, ...]:
    """Parse ``8x6,16x12`` (patch) or ``chord=0.01,chord=0.001`` (trimmed)."""
    rungs = []
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if token.startswith("chord="):
            rungs.append(Rung(chord_tolerance=float(token.split("=", 1)[1])))
        else:
            u, _, v = token.partition("x")
            rungs.append(Rung(u_samples=int(u), v_samples=int(v)))
    if not rungs:
        raise ValueError(f"no rungs in {spec!r}")
    return tuple(rungs)


# ---------------------------------------------------------------------------
# measurement helpers
# ---------------------------------------------------------------------------
def _fibonacci_directions(count: int) -> list[tuple[float, float, float]]:
    """Deterministic, roughly uniform ray directions (same set as the generator ladder)."""
    golden = math.pi * (3.0 - math.sqrt(5.0))
    directions = []
    for index in range(count):
        z = 1.0 - 2.0 * (index + 0.5) / count
        radius = math.sqrt(max(0.0, 1.0 - z * z))
        theta = golden * index
        directions.append((math.cos(theta) * radius, math.sin(theta) * radius, z))
    return directions


def _finite(value: Any) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def mesh_extent(vertices: Sequence[Sequence[float]]) -> float:
    """Largest absolute coordinate of the model: the scale the weld uses."""
    extent = 0.0
    for vertex in vertices:
        for coordinate in vertex:
            extent = max(extent, abs(float(coordinate)))
    return extent


def _open_boundary_count(mesh: Any) -> int:
    from .mesh import boundary_edges

    return len(boundary_edges(mesh))


def _watertight(mesh: Any) -> bool:
    from .mesh import is_watertight

    return bool(is_watertight(mesh))


def read_source(path: str | Path, route: str) -> tuple[Any, float, dict[str, Any]]:
    """Read the CAD file for ``route``; returns ``(model, seconds, census)``."""
    started = time.perf_counter()
    if route == ROUTE_PATCH:
        from .geomquality.io import read_iges, read_step

        suffix = Path(path).suffix.lower()
        reader = read_step if suffix in (".step", ".stp") else read_iges
        model = reader(path)
        census = {
            "faces": len(model.faces),
            "edges": len(model.edges),
            "surfaces": len(model.faces),
        }
    elif route == ROUTE_TRIMMED:
        from .geomquality.trim_mesh import read_iges_trimmed

        model = read_iges_trimmed(path)
        census = {
            "faces": len(model.faces),
            "metres_per_unit": model.metres_per_unit,
        }
    else:
        raise ValueError(f"unknown route {route!r}; expected one of {ROUTES}")
    return model, time.perf_counter() - started, census


def tessellate_rung(model: Any, route: str, rung: Rung, *,
                    chunk_faces: int | None = None,
                    on_chunk: Any = None) -> tuple[Any, float]:
    """Tessellate one rung; returns ``(mesh, seconds)``."""
    started = time.perf_counter()
    if route == ROUTE_PATCH:
        mesh = model.tessellate(rung.u_samples, rung.v_samples,
                                chunk_faces=chunk_faces, on_chunk=on_chunk)
    else:
        from .geomquality.trim_mesh import tessellate_trimmed_model

        assert rung.chord_tolerance is not None
        tessellation = tessellate_trimmed_model(
            model, scale=model.metres_per_unit,
            chord_tolerance=rung.chord_tolerance, refine_levels=6)
        from .geometry import Mesh

        mesh = Mesh(vertices=tuple(tessellation.vertices), faces=tuple(tessellation.faces))
    return mesh, time.perf_counter() - started


def _result(subtask: str, rung: Rung, mesh: Any, seconds: float, *,
            route: str, note: str = "", **extra: Any) -> SubTaskResult:
    triangles = len(mesh.faces) if mesh is not None else 0
    vertices = len(mesh.vertices) if mesh is not None else 0
    return SubTaskResult(
        subtask=subtask,
        nx=rung.u_samples,
        nz=rung.v_samples,
        triangles=triangles,
        vertices=vertices,
        grid_cells=0,
        seconds=round(seconds, 6),
        status=MEASURED,
        throughput=round(triangles / seconds, 1) if seconds > 0 and triangles else None,
        note=note,
        peak_rss_mb=peak_rss_mb(),
        rss_mb=live_rss_mb(),
        mesh_mb=mesh_footprint_mb(mesh) if mesh is not None else None,
        **extra,
    )


# ---------------------------------------------------------------------------
# one rung
# ---------------------------------------------------------------------------
def measure_rung(path: str | Path, route: str, rung: Rung, *, model: Any = None,
                 digest: str = "", queries: int = DEFAULT_QUERIES,
                 chunk_faces: int | None = DEFAULT_CHUNK_FACES,
                 max_triangles: int = DEFAULT_MAX_TRIANGLES) -> tuple[dict[str, Any], Any]:
    """Measure every sub-task on one rung; returns ``(record, mesh)``.

    The record is a plain dict (it is what lands in the artifact); the mesh is
    returned so a caller can skip re-tessellating for the next sub-task.
    """
    parameter = rung.parameter(route)
    if model is None:
        model, parse_seconds, census = read_source(path, route)
    else:
        parse_seconds, census = None, {}

    face_count = len(getattr(model, "faces", ()) or ())
    expected = rung.expected_triangles(face_count) if route == ROUTE_PATCH else 0
    if expected and expected > max_triangles:
        record = {
            "subtask": "cad_tessellate", "route": route, "parameter": parameter,
            "status": NOT_RUN, "seconds": None, "triangles": 0, "vertices": 0,
            "note": (f"not_run: expected {expected:,} triangles exceeds "
                     f"max_triangles={max_triangles:,}; raise --max-triangles to run it"),
        }
        return record, None

    mesh, tess_seconds = tessellate_rung(model, route, rung)

    record: dict[str, Any] = {
        "subtask": "cad_tessellate",
        "route": route,
        "parameter": parameter,
        "status": MEASURED,
        "seconds": round(tess_seconds, 6),
        "triangles": len(mesh.faces),
        "vertices": len(mesh.vertices),
        "throughput": round(len(mesh.faces) / tess_seconds, 1) if tess_seconds > 0 else None,
        "parse_seconds": round(parse_seconds, 6) if parse_seconds is not None else None,
        "census": census,
        "mesh_mb": mesh_footprint_mb(mesh),
        "rss_mb": live_rss_mb(),
        "peak_rss_mb": peak_rss_mb(),
        "note": "",
    }

    # --- chunked tessellation: same geometry, bounded sampling peak ---------
    if chunk_faces and route == ROUTE_PATCH:
        chunk_rss: list[float] = []

        def _on_chunk(done: int, total: int) -> None:
            value = live_rss_mb()
            if value is not None:
                chunk_rss.append(round(value, 1))

        chunked, chunk_seconds = tessellate_rung(model, route, rung,
                                                 chunk_faces=chunk_faces,
                                                 on_chunk=_on_chunk)
        identical = (tuple(chunked.vertices) == tuple(mesh.vertices)
                     and tuple(chunked.faces) == tuple(mesh.faces))
        record["chunked"] = {
            "chunk_faces": chunk_faces,
            "seconds": round(chunk_seconds, 6),
            "identical": bool(identical),
            "rss_per_chunk_mb": chunk_rss,
            "rss_peak_of_chunks_mb": max(chunk_rss) if chunk_rss else None,
            "rss_first_chunk_mb": chunk_rss[0] if chunk_rss else None,
            "speed_ratio": (round(chunk_seconds / tess_seconds, 3)
                            if tess_seconds > 0 else None),
        }

    # --- tolerance welding --------------------------------------------------
    from .mesh import weld

    extent = mesh_extent(mesh.vertices)
    tolerance = DEFAULT_WELD_RELATIVE_TOLERANCE * max(extent, 1e-12)
    boundary_before = _open_boundary_count(mesh)
    started = time.perf_counter()
    welded = weld(mesh, tolerance)
    weld_seconds = time.perf_counter() - started
    boundary_after = _open_boundary_count(welded)
    record["weld"] = {
        "tolerance": tolerance,
        "seconds": round(weld_seconds, 6),
        "vertices_before": len(mesh.vertices),
        "vertices_after": len(welded.vertices),
        "vertices_merged": len(mesh.vertices) - len(welded.vertices),
        "merged_fraction": (round((len(mesh.vertices) - len(welded.vertices))
                                  / len(mesh.vertices), 6) if mesh.vertices else None),
        "boundary_edges_before": boundary_before,
        "boundary_edges_after": boundary_after,
        "watertight_before": _watertight(mesh),
        "watertight_after": _watertight(welded),
        "triangles": len(welded.faces),
    }

    # --- volume -------------------------------------------------------------
    from .mesh import total_signed_volume

    started = time.perf_counter()
    volume = total_signed_volume(welded)
    volume_seconds = time.perf_counter() - started
    from .mesh import surface_area

    record["volume"] = {
        "seconds": round(volume_seconds, 6),
        "signed_volume": volume,
        "abs_volume": abs(volume),
        "surface_area": surface_area(welded),
        "watertight": _watertight(welded),
    }

    # --- accelerated point location ----------------------------------------
    from .spatial_index import TriangleGrid

    started = time.perf_counter()
    grid = TriangleGrid(welded)
    build_seconds = time.perf_counter() - started
    directions = _fibonacci_directions(queries)
    centroid = _centroid(welded.vertices)
    started = time.perf_counter()
    hits = sum(1 for direction in directions
               if grid.cast_ray(centroid, direction) is not None)
    query_seconds = time.perf_counter() - started
    stats = grid.stats()
    record["locate"] = {
        "queries": queries,
        "seconds": round(query_seconds, 6),
        "queries_per_second": round(queries / query_seconds, 1) if query_seconds > 0 else None,
        "build_seconds": round(build_seconds, 6),
        "index_cells": int(stats["cells"]),
        "index_entries": int(stats["entries"]),
        "index_bytes": int(stats["index_bytes"]),
        "hits": hits,
    }

    record["rss_mb_end"] = live_rss_mb()
    record["peak_rss_mb_end"] = peak_rss_mb()
    return record, welded


def _centroid(vertices: Sequence[Sequence[float]]) -> tuple[float, float, float]:
    count = len(vertices) or 1
    sx = sum(float(v[0]) for v in vertices) / count
    sy = sum(float(v[1]) for v in vertices) / count
    sz = sum(float(v[2]) for v in vertices) / count
    return (sx, sy, sz)


# ---------------------------------------------------------------------------
# regression gate (2x = RED)
# ---------------------------------------------------------------------------
def regression_gates(records: Iterable[dict[str, Any]], baseline: dict[str, Any] | None,
                     digest: str, factor: float = DEFAULT_GATE_FACTOR) -> tuple[GateVerdict, ...]:
    """Wall-clock gate: a rung more than ``factor`` slower than its baseline is RED.

    A rung with no baseline entry is ``not_run`` with the reason spelled out -
    an absent baseline is not evidence of health.  Chunked and weld times are
    gated on their own keys, because a regression in the welding stage would be
    invisible in the tessellation time.
    """
    stored: dict[str, float] = {}
    if baseline:
        for key, value in (baseline.get("measurements") or {}).items():
            if _finite(value):
                stored[key] = float(value)

    gates: list[GateVerdict] = []
    for record in records:
        route = record["route"]
        parameter = record["parameter"]
        pairs = [
            ("cad_tessellate", record.get("seconds")),
            ("cad_weld", (record.get("weld") or {}).get("seconds")),
            ("cad_volume", (record.get("volume") or {}).get("seconds")),
            ("cad_locate", (record.get("locate") or {}).get("seconds")),
        ]
        if record.get("chunked"):
            pairs.append(("cad_tessellate_chunked", record["chunked"].get("seconds")))
        for subtask, seconds in pairs:
            name = f"regression_{subtask}[{route} {parameter}]"
            key = measurement_key(subtask, digest, route, parameter)
            if seconds is None:
                gates.append(GateVerdict(
                    name=name, passed=True, status=NOT_RUN,
                    detail="not_run: no measurement for this rung (see the rung note)"))
                continue
            if key not in stored:
                gates.append(GateVerdict(
                    name=name, passed=True, status=NOT_RUN,
                    detail=f"not_run: no baseline entry {key!r}; run --save-baseline first"))
                continue
            base = stored[key]
            ratio = float(seconds) / base if base > 0 else math.inf
            passed = ratio <= factor
            gates.append(GateVerdict(
                name=name, passed=passed, status=MEASURED,
                detail=(f"{ratio:.2f}x of baseline {base:.4f}s "
                        f"(limit {factor:.1f}x) -> {'ok' if passed else 'RED'}")))
    return tuple(gates)


def chunk_equivalence_gates(records: Iterable[dict[str, Any]]) -> tuple[GateVerdict, ...]:
    """Chunking must not change the mesh; a mismatch is a correctness RED."""
    gates = []
    for record in records:
        chunked = record.get("chunked")
        if not chunked:
            continue
        name = f"chunk_equivalence[{record['route']} {record['parameter']}]"
        same = bool(chunked.get("identical"))
        gates.append(GateVerdict(
            name=name, passed=same, status=MEASURED,
            detail=("chunked tessellation is bit-identical to the one-shot mesh"
                    if same else "RED: chunked tessellation changed the mesh")))
    return tuple(gates)


# ---------------------------------------------------------------------------
# ladder
# ---------------------------------------------------------------------------
def _environment() -> dict[str, Any]:
    try:
        load1, load5, load15 = os.getloadavg()
    except OSError:  # pragma: no cover - platform without loadavg
        load1 = load5 = load15 = None
    from . import __version__ as tensorcad_version

    return {
        "tensorcad": tensorcad_version,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "loadavg": {"load1": load1, "load5": load5, "load15": load15},
    }


def run_cad_ladder(path: str | Path, rungs: Sequence[Rung], *, route: str = ROUTE_PATCH,
                   queries: int = DEFAULT_QUERIES,
                   chunk_faces: int | None = DEFAULT_CHUNK_FACES,
                   max_triangles: int = DEFAULT_MAX_TRIANGLES,
                   baseline: dict[str, Any] | None = None,
                   gate_factor: float = DEFAULT_GATE_FACTOR) -> dict[str, Any]:
    """Run one rung of every sub-task per rung and gate the result."""
    started = time.time()
    digest = source_digest(path)
    model, parse_seconds, census = read_source(path, route)
    records: list[dict[str, Any]] = []
    meshes: list[Any] = []
    for rung in rungs:
        record, mesh = measure_rung(path, route, rung, model=model, digest=digest,
                                    queries=queries, chunk_faces=chunk_faces,
                                    max_triangles=max_triangles)
        record["parse_seconds"] = round(parse_seconds, 6)
        records.append(record)
        if mesh is not None:
            meshes.append(mesh)

    gates = list(chunk_equivalence_gates(records))
    gates.extend(regression_gates(records, baseline, digest, gate_factor))

    # Scaling exponent: reuse the generated ladder's implementation by mapping
    # each rung onto the sub-task results it already understands.
    scaling_results = []
    for record in records:
        if record["status"] != MEASURED:
            continue
        rung = _rung_from_record(record, route)
        triangles = record["triangles"]
        for subtask, seconds in (
            ("cad_tessellate", record.get("seconds")),
            ("cad_weld", (record.get("weld") or {}).get("seconds")),
            ("cad_locate", (record.get("locate") or {}).get("seconds")),
        ):
            if seconds is None:
                continue
            scaling_results.append(SubTaskResult(
                subtask=subtask, nx=rung.u_samples, nz=rung.v_samples,
                triangles=triangles, vertices=record["vertices"], grid_cells=0,
                seconds=float(seconds), status=MEASURED))
    gates.extend(scaling_gates(
        scaling_results,
        limits={"cad_tessellate": 1.5, "cad_weld": 1.5, "cad_locate": 1.5}))

    payload = {
        "schema": SCHEMA,
        "route": route,
        "source": {
            "path": str(path),
            "sha256": digest,
            "bytes": Path(path).stat().st_size,
            "parse_seconds": round(parse_seconds, 6),
            "census": census,
        },
        "environment": _environment(),
        "gate_factor": gate_factor,
        "queries": queries,
        "chunk_faces": chunk_faces,
        "max_triangles": max_triangles,
        "wall_seconds": round(time.time() - started, 3),
        "records": records,
        "measurements": _measurement_table(records, digest),
        "gates": [dict(name=g.name, passed=g.passed, status=g.status, detail=g.detail)
                  for g in gates],
        "cluster_thresholds": cluster_thresholds(records),
    }
    payload["passed"] = all(g.passed for g in gates)
    payload["gates_checked"] = sum(1 for g in gates if g.status != NOT_RUN)
    payload["gates_not_run"] = sum(1 for g in gates if g.status == NOT_RUN)
    payload["verdict"] = _verdict(gates)
    return payload


def _rung_from_record(record: dict[str, Any], route: str) -> Rung:
    if route == ROUTE_TRIMMED:
        tolerance = float(str(record["parameter"]).split("=", 1)[1])
        return Rung(chord_tolerance=tolerance)
    u, _, v = str(record["parameter"]).split("=", 1)[1].partition("x")
    return Rung(u_samples=int(u), v_samples=int(v))


def _measurement_table(records: Iterable[dict[str, Any]], digest: str) -> dict[str, float]:
    """``key -> seconds``: exactly what the regression gate reads back."""
    table: dict[str, float] = {}
    for record in records:
        if record["status"] != MEASURED:
            continue
        pairs = [
            ("cad_tessellate", record.get("seconds")),
            ("cad_weld", (record.get("weld") or {}).get("seconds")),
            ("cad_volume", (record.get("volume") or {}).get("seconds")),
            ("cad_locate", (record.get("locate") or {}).get("seconds")),
        ]
        if record.get("chunked"):
            pairs.append(("cad_tessellate_chunked", record["chunked"].get("seconds")))
        for subtask, seconds in pairs:
            if seconds is None:
                continue
            table[measurement_key(subtask, digest, record["route"],
                                  record["parameter"])] = float(seconds)
    return table


def cluster_thresholds(records: Iterable[dict[str, Any]],
                       budget_seconds: float = 60.0) -> dict[str, Any]:
    """At what triangle count does each sub-task leave a ``budget_seconds`` budget?

    Fit a log-log slope through the measured rungs of each sub-task and invert
    it: the result is "the size at which this stage stops fitting in the
    budget", the reproducible form of the anchor's "record at what scale the
    cluster is needed".  Fewer than two positive points -> ``unavailable`` with
    the reason, never a guess.
    """
    grouped: dict[str, list[tuple[float, float]]] = {}
    for record in records:
        if record["status"] != MEASURED or not record["triangles"]:
            continue
        size = float(record["triangles"])
        for subtask, seconds in (
            ("cad_tessellate", record.get("seconds")),
            ("cad_weld", (record.get("weld") or {}).get("seconds")),
            ("cad_volume", (record.get("volume") or {}).get("seconds")),
            ("cad_locate", (record.get("locate") or {}).get("seconds")),
        ):
            if seconds is None or float(seconds) <= 0.0:
                continue
            grouped.setdefault(subtask, []).append((size, float(seconds)))

    out: dict[str, Any] = {}
    for subtask, points in sorted(grouped.items()):
        points.sort()
        if len(points) < 2:
            out[subtask] = {"status": NOT_RUN,
                            "reason": "fewer than 2 measured sizes: no slope to fit"}
            continue
        (x0, y0), (x1, y1) = points[0], points[-1]
        if x1 <= x0 or y0 <= 0 or y1 <= 0:
            out[subtask] = {"status": NOT_RUN, "reason": "degenerate measurements"}
            continue
        exponent = (math.log(y1) - math.log(y0)) / (math.log(x1) - math.log(x0))
        rate = y1 / x1 ** exponent
        if rate <= 0:
            out[subtask] = {"status": NOT_RUN, "reason": "non-positive rate"}
            continue
        size_at_budget = (budget_seconds / rate) ** (1.0 / exponent)
        out[subtask] = {
            "status": MEASURED,
            "exponent": round(exponent, 3),
            "budget_seconds": budget_seconds,
            "size_at_budget": round(size_at_budget, 1),
            "measured_max_size": int(x1),
            "measured_max_seconds": round(y1, 4),
            "exceeds_max_rung": bool(size_at_budget > x1),
        }
    return out


def _verdict(gates: Sequence[GateVerdict]) -> str:
    """``fail`` / ``not_run`` / ``pass`` for the ladder as a whole.

    A correctness gate (chunk equivalence) passing is **not** a performance
    verdict: with no wall-clock baseline every regression gate is ``not_run``,
    and calling that a pass would report health that was never measured.  So a
    pass needs at least one evaluated ``regression_*`` gate.
    """
    if any(not g.passed for g in gates):
        return "fail"
    measured = [g for g in gates if g.status != NOT_RUN]
    if not any(g.name.startswith("regression_") for g in measured):
        return "not_run"
    return "pass"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _print_table(payload: dict[str, Any]) -> None:
    print(f"source: {payload['source']['path']}  sha256={payload['source']['sha256'][:16]}  "
          f"route={payload['route']}")
    census = payload["source"]["census"]
    print(f"census: {census}")
    header = (f"{'chord':>8} {'uv':>10} {'triangles':>12} {'vertices':>12} "
              f"{'tess s':>8} {'weld s':>8} {'weld v_merged':>14} {'bnd->':>7} {'locate q/s':>11}")
    print(header)
    for record in payload["records"]:
        if record["status"] != MEASURED:
            print(f"{record['parameter']:>19} {'not_run':>12}  {record['note']}")
            continue
        weld = record.get("weld") or {}
        locate = record.get("locate") or {}
        print(f"{record['parameter']:>19} {record['triangles']:>12} {record['vertices']:>12} "
              f"{record['seconds']:>8.3f} {weld.get('seconds', float('nan')):>8.3f} "
              f"{weld.get('vertices_merged', 0):>14} "
              f"{weld.get('boundary_edges_before', 0)}->{weld.get('boundary_edges_after', 0)} "
              f"{locate.get('queries_per_second', float('nan')):>10.0f}")
    print()
    for entry in payload["gates"]:
        state = "ok" if entry["passed"] else "RED"
        if entry["status"] == NOT_RUN:
            state = "not_run"
        print(f"  [{state:>7}] {entry['name']}: {entry['detail']}")
    print(f"\nverdict: {payload['verdict']}  "
          f"(gates checked={payload['gates_checked']} not_run={payload['gates_not_run']})")
    if payload["verdict"] == "not_run":
        print("This is NOT a pass: no gate could be evaluated.")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m tensorcad.perf_cad_source",
        description="Performance ladder on a real IGES/STEP source (GAP-CAD-D4-PERF-SCALE).")
    parser.add_argument("--source", default=DEFAULT_PATCH_SOURCE)
    parser.add_argument("--route", choices=ROUTES, default=ROUTE_PATCH)
    parser.add_argument("--rungs", default=None,
                        help="patch: 8x6,16x12 ; trimmed: chord=0.01,chord=0.001")
    parser.add_argument("--quick", action="store_true",
                        help="two small rungs (smoke test, seconds not minutes)")
    parser.add_argument("--queries", type=int, default=DEFAULT_QUERIES)
    parser.add_argument("--chunk-faces", type=int, default=DEFAULT_CHUNK_FACES)
    parser.add_argument("--no-chunk", action="store_true", help="skip the chunked sub-task")
    parser.add_argument("--max-triangles", type=int, default=DEFAULT_MAX_TRIANGLES)
    parser.add_argument("--baseline", default=None, help="baseline JSON to gate against")
    parser.add_argument("--save-baseline", default=None, help="write the baseline JSON here")
    parser.add_argument("--regression-factor", type=float, default=DEFAULT_GATE_FACTOR)
    parser.add_argument("--out", default=None, help="write the full artifact JSON here")
    args = parser.parse_args(argv)

    if args.quick:
        rungs = ((Rung(u_samples=4, v_samples=3), Rung(u_samples=8, v_samples=6))
                 if args.route == ROUTE_PATCH
                 else (Rung(chord_tolerance=5e-2), Rung(chord_tolerance=2e-2)))
    elif args.rungs:
        rungs = parse_rungs(args.rungs)
    else:
        rungs = (tuple(Rung(u_samples=u, v_samples=v) for u, v in DEFAULT_SAMPLE_RUNGS)
                 if args.route == ROUTE_PATCH
                 else tuple(Rung(chord_tolerance=t) for t in DEFAULT_CHORD_TOLERANCES))

    baseline = None
    if args.baseline:
        baseline = json.loads(Path(args.baseline).read_text())

    payload = run_cad_ladder(
        args.source, rungs, route=args.route, queries=args.queries,
        chunk_faces=None if args.no_chunk else args.chunk_faces,
        max_triangles=args.max_triangles, baseline=baseline,
        gate_factor=args.regression_factor)

    _print_table(payload)

    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=2) + "\n")
        print(f"\nartifact: {out}")
    if args.save_baseline:
        target = Path(args.save_baseline)
        target.parent.mkdir(parents=True, exist_ok=True)
        previous = {}
        if target.exists():
            previous = (json.loads(target.read_text()).get("measurements") or {})
        previous.update(payload["measurements"])
        target.write_text(json.dumps({
            "schema": SCHEMA,
            "generated_from": payload["source"]["path"],
            "source_sha256": payload["source"]["sha256"],
            "environment": payload["environment"],
            "baseline_policy": "min of the ladder's own rungs; one key per "
                              "(subtask, source digest, route, parameter)",
            "measurements": previous,
        }, indent=2, sort_keys=True) + "\n")
        print(f"baseline: {target} ({len(previous)} keys)")
    return 0 if payload["verdict"] in ("pass", "not_run") else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())