#!/usr/bin/env python3
"""Benchmark the OpenFOAM polyMesh route of ``tensorcad.meshquality`` on real ship geometry.

The numbers, and how each can fail:

1. **Closed forms** -- a uniformly sheared Cartesian grid must give ``atan(shear)``
   non-orthogonality and *zero* skewness; two stacked cells of different extent must give
   45 deg and ``0.5/sqrt(2)``. Analytic values, tolerance 1e-12.
2. **Independent route** -- :func:`tensorcad.meshquality.topology.non_orthogonality` on the
   same cells. Identical on affine meshes (parallelepipeds, where the vertex average is the
   centroid); on the hull mesh the two differ, and the report says by how much and why.
3. **Round trip** -- every mesh goes through ``write_polymesh`` -> ``read_polymesh`` and the
   measurements must be bit-identical to the in-memory ones.
4. **Real geometry** -- the repository's hull volume mesh, measured through the solver's own
   container: non-orthogonality, skewness, and the collapsed keel/stem cells the builder
   documents.

Writes ``artifacts/mesh_quality_polymesh/results.json`` and ``plot.png``; exits non-zero if
an analytic or round-trip criterion stops holding.
"""
from __future__ import annotations

import json
import math
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tensorcad.hull import HullParameters  # noqa: E402
from tensorcad.meshquality import topology  # noqa: E402
from tensorcad.meshquality.polymesh import (  # noqa: E402
    DEFAULT_POLYMESH_LIMITS,
    check_polymesh,
    polymesh_from_hexes,
    read_polymesh,
    write_polymesh,
)
from tensorcad.volume_mesh import volume_mesh_from_hull  # noqa: E402

ARTIFACT_DIR = REPO_ROOT / "artifacts" / "mesh_quality_polymesh"
FAILURES: list[str] = []


def box_grid(nx: int, ny: int, nz: int):
    points, index = [], {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                index[(i, j, k)] = len(points)
                points.append((float(i), float(j), float(k)))
    hexes = [
        [
            index[(i, j, k)], index[(i + 1, j, k)], index[(i + 1, j + 1, k)], index[(i, j + 1, k)],
            index[(i, j, k + 1)], index[(i + 1, j, k + 1)], index[(i + 1, j + 1, k + 1)],
            index[(i, j + 1, k + 1)],
        ]
        for i in range(nx)
        for j in range(ny)
        for k in range(nz)
    ]
    return points, hexes


def measure(name: str, points, hexes) -> dict:
    mesh = polymesh_from_hexes(points, hexes)
    quality = check_polymesh(mesh)
    with tempfile.TemporaryDirectory() as directory:
        write_polymesh(mesh, directory)
        again = check_polymesh(read_polymesh(directory))
    roundtrip = (
        again.non_orthogonality == quality.non_orthogonality
        and again.skewness == quality.skewness
        and again.volumes == quality.volumes
    )
    if not roundtrip:
        FAILURES.append(f"{name}: polyMesh round trip changed the measurements")
    ortho = [a for a in quality.non_orthogonality if a == a]
    skew = [s for s in quality.skewness if s == s]
    hex_route = topology.non_orthogonality(points, hexes, "hex").maximum
    record = {
        "mesh": name,
        "cells": mesh.n_cells,
        "internal_faces": mesh.n_internal_faces,
        "boundary_faces": mesh.n_boundary_faces,
        "max_non_orthogonality": max(ortho) if ortho else None,
        "mean_non_orthogonality": sum(ortho) / len(ortho) if ortho else None,
        "max_skewness": max(skew) if skew else None,
        "min_cell_volume": min(quality.volumes),
        "total_volume": math.fsum(quality.volumes),
        "non_positive_cells": sum(1 for v in quality.volumes if v <= 0.0),
        "degenerate_faces": sum(1 for a in quality.non_orthogonality if a != a),
        "hex_route_max_non_orthogonality": hex_route,
        "route_difference_deg": (
            (max(ortho) - hex_route) if ortho and hex_route == hex_route else None
        ),
        "roundtrip_identical": roundtrip,
        "gate_passed": quality.gate().passed,
        "gate_failures": [c.name for c in quality.gate().failures],
        "non_orthogonality": ortho,
        "skewness": skew,
    }
    print(
        f"{name:22s} cells={record['cells']:5d} internal={record['internal_faces']:5d} "
        f"maxNonOrth={fmt(record['max_non_orthogonality'])} maxSkew={fmt(record['max_skewness'])} "
        f"minV={record['min_cell_volume']:.3e} gate={'PASS' if record['gate_passed'] else 'FAIL'}"
    )
    return record


def fmt(value):
    return "   (nan)" if value is None else f"{value:8.4f}"


def main() -> int:
    results: dict[str, object] = {"limits": DEFAULT_POLYMESH_LIMITS, "analytic": [], "meshes": []}

    # 1. closed forms -------------------------------------------------------------------
    points, hexes = box_grid(3, 3, 3)
    shear = 0.5
    sheared = measure("sheared_box", [(x + shear * y, y, z) for (x, y, z) in points], hexes)
    expected_angle = math.degrees(math.atan(shear))
    analytic = {
        "case": "uniform_shear_0.5",
        "non_orthogonality": sheared["max_non_orthogonality"],
        "expected": expected_angle,
        "difference": abs(sheared["max_non_orthogonality"] - expected_angle),
        "skewness": sheared["max_skewness"],
        "expected_skewness": 0.0,
    }
    if analytic["difference"] > 1e-12 or abs(analytic["skewness"]) > 1e-12:
        FAILURES.append(f"uniform shear: {analytic}")
    results["analytic"].append(analytic)

    stacked_points = [
        (0, 0, -1), (1, 0, -1), (1, 1, -1), (0, 1, -1),
        (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
        (2, 0, 1), (3, 0, 1), (3, 1, 1), (2, 1, 1),
    ]
    stacked = measure("stacked_shear", stacked_points, [list(range(8)), [4, 5, 6, 7, 8, 9, 10, 11]])
    expected_skew = 0.5 / math.sqrt(2.0)
    stacked_case = {
        "case": "stacked_cells_shifted_by_2",
        "non_orthogonality": stacked["max_non_orthogonality"],
        "expected": 45.0,
        "skewness": stacked["max_skewness"],
        "expected_skewness": expected_skew,
        "skewness_difference": abs(stacked["max_skewness"] - expected_skew),
    }
    if abs(stacked_case["non_orthogonality"] - 45.0) > 1e-12 or stacked_case["skewness_difference"] > 1e-12:
        FAILURES.append(f"stacked cells: {stacked_case}")
    results["analytic"].append(stacked_case)

    box_points, box_hexes = box_grid(4, 3, 2)
    box = measure("cartesian_box_4x3x2", box_points, box_hexes)
    if box["max_non_orthogonality"] > 1e-12 or box["max_skewness"] > 1e-12:
        FAILURES.append(f"cartesian box is not orthogonal: {box['max_non_orthogonality']}")

    # 2. the repository's own ship geometry ---------------------------------------------
    meshes = [box, sheared, stacked]
    for layers, nx, nz in ((1, 17, 9), (3, 33, 17)):
        hull = volume_mesh_from_hull(
            HullParameters(nx=nx, nz=nz), first_layer=0.05, layers=layers
        )
        meshes.append(measure(f"hull_nx{nx}_nz{nz}_L{layers}", hull.vertices, hull.hexes))
    results["meshes"] = [
        {key: value for key, value in record.items() if key not in ("non_orthogonality", "skewness")}
        for record in meshes
    ]

    # 3. plot ---------------------------------------------------------------------------
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    hull_records = [record for record in meshes if record["mesh"].startswith("hull")]
    for record in hull_records:
        axes[0].hist(record["non_orthogonality"], bins=40, alpha=0.6, label=record["mesh"])
    axes[0].axvline(DEFAULT_POLYMESH_LIMITS["max_non_orthogonality"], color="r", ls="--",
                    label=f"checkMesh limit {DEFAULT_POLYMESH_LIMITS['max_non_orthogonality']:g} deg")
    axes[0].set_xlabel("internal-face non-orthogonality [deg]")
    axes[0].set_ylabel("faces")
    axes[0].set_title("Hull volume mesh through polyMesh")
    axes[0].legend(fontsize=7)

    for record in hull_records:
        axes[1].hist(record["skewness"], bins=40, alpha=0.6, label=record["mesh"])
    axes[1].axvline(DEFAULT_POLYMESH_LIMITS["max_skewness"], color="r", ls="--",
                    label=f"checkMesh limit {DEFAULT_POLYMESH_LIMITS['max_skewness']:g}")
    axes[1].set_xlabel("face skewness [-]")
    axes[1].set_ylabel("faces")
    axes[1].set_title("Skewness (0 on a parallelepiped)")
    axes[1].legend(fontsize=7)

    labels = [record["mesh"] for record in meshes]
    poly = [record["max_non_orthogonality"] or 0.0 for record in meshes]
    hexroute = [record["hex_route_max_non_orthogonality"] or 0.0 for record in meshes]
    positions = range(len(labels))
    axes[2].bar([p - 0.2 for p in positions], poly, width=0.4, label="polyMesh route (OpenFOAM centres)")
    axes[2].bar([p + 0.2 for p in positions], hexroute, width=0.4, label="hex route (vertex-average centres)")
    axes[2].axhline(math.degrees(math.atan(0.5)), color="k", ls=":", label="atan(0.5) analytic")
    axes[2].set_xticks(list(positions))
    axes[2].set_xticklabels(labels, rotation=30, ha="right", fontsize=7)
    axes[2].set_ylabel("max non-orthogonality [deg]")
    axes[2].set_title("Two routes; identical on affine meshes")
    axes[2].legend(fontsize=7)
    figure.tight_layout()
    figure.savefig(ARTIFACT_DIR / "plot.png", dpi=140)
    plt.close(figure)

    results["failures"] = FAILURES
    (ARTIFACT_DIR / "results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"\nartifacts: {ARTIFACT_DIR}")
    if FAILURES:
        for failure in FAILURES:
            print(f"FAIL {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())