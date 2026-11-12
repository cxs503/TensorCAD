#!/usr/bin/env python3
"""Benchmark and self-check for the SDF voxel-grid export (``to_voxel_grid``).

What it measures, in order:

1. **Accuracy vs resolution** for a sphere of known volume: the partial-fill volume
   (second-order cut rule) against hard cell counting (what a boolean occupancy grid
   gives you), plus the band-count area estimate.
2. **Anisotropic cells** on an ellipsoid, where one cell per axis is what keeps a
   4:2:1 body from spending its grid on empty water.
3. **The published hull**: the SUBOFF bare hull against the displacement volume of
   the DARPA offset table (0.699211 m^3), so the voxel path answers to the same
   external reference the mesh path does.
4. **Solver handoff**: the record a consumer reads, including its SHA-256.

Exits non-zero if any claim stops holding, so CI and the docs cannot drift apart.
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_REPO_ROOT = _Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in _sys.path:  # allow running from a checkout without installing
    _sys.path.insert(0, str(_REPO_ROOT))

import argparse
import json
import math
import time
from pathlib import Path

from tensorcad import EllipsoidField, SphereField, SuboffHullField, to_voxel_grid

SPHERE_RADIUS = 0.5
SPHERE_VOLUME = 4.0 / 3.0 * math.pi * SPHERE_RADIUS**3
SUBOFF_VOLUME = 0.699211
SUBOFF_WETTED_AREA = 5.988264
ELLIPSOID_SEMI_AXES = (0.5, 0.25, 0.125)
ELLIPSOID_VOLUME = 4.0 / 3.0 * math.pi * ELLIPSOID_SEMI_AXES[0] * ELLIPSOID_SEMI_AXES[1] * ELLIPSOID_SEMI_AXES[2]

BOX = ((-0.6, -0.6, -0.6), (0.6, 0.6, 0.6))


def _fit_order(points: list[tuple[int, float]]) -> float:
    """Least-squares slope of log|error| against log(cells): the observed order."""
    xs = [math.log(cells) for cells, _ in points]
    ys = [math.log(abs(error)) for _, error in points]
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    numerator = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys, strict=True))
    denominator = sum((x - mean_x) ** 2 for x in xs)
    # |error| shrinks with resolution, so the fit of log|e| vs log(h) is negative; the
    # convergence *order* is its magnitude.
    return -numerator / denominator


def sphere_table(resolutions=(8, 16, 32, 64)) -> dict[str, object]:
    print("## Sphere, r = 0.5 m, box +-0.6 m (non-coincident with the surface)\n")
    print(f"exact volume {SPHERE_VOLUME:.9f} m^3\n")
    print("| cells | cells total | partial-fill volume [m^3] | err % | hard-count volume [m^3] | err % |")
    print("|---|---|---|---|---|---|")
    partial_points: list[tuple[int, float]] = []
    hard_points: list[tuple[int, float]] = []
    rows = []
    for cells in resolutions:
        grid = to_voxel_grid(SphereField(radius=SPHERE_RADIUS), cells, bounds=BOX)
        partial = grid.volume()
        hard = grid.occupied_volume()
        partial_error = partial / SPHERE_VOLUME - 1.0
        hard_error = hard / SPHERE_VOLUME - 1.0
        partial_points.append((cells, partial_error))
        hard_points.append((cells, hard_error))
        rows.append(
            {
                "cells": cells,
                "volume_partial": partial,
                "error_partial": partial_error,
                "volume_hard": hard,
                "error_hard": hard_error,
            }
        )
        print(
            f"| {cells} | {grid.cells:,} | {partial:.9f} | {100 * partial_error:+.4f} | "
            f"{hard:.9f} | {100 * hard_error:+.4f} |"
        )
    order_partial = _fit_order(partial_points)
    order_hard = _fit_order(hard_points)
    print(f"\nObserved convergence order (least squares over the table): "
          f"partial-fill {order_partial:.2f}, hard count {order_hard:.2f}.\n")
    assert abs(rows[-1]["error_partial"]) < 1e-3, rows[-1]
    assert abs(rows[-1]["error_partial"]) < abs(rows[-1]["error_hard"]) / 2.0, rows[-1]
    assert order_partial > 1.0, order_partial
    return {"rows": rows, "order_partial": order_partial, "order_hard": order_hard}


def ellipsoid_table() -> dict[str, object]:
    print("## Ellipsoid, semi-axes (0.5, 0.25, 0.125) m: cells per axis\n")
    print(f"exact volume {ELLIPSOID_VOLUME:.9f} m^3\n")
    print("| cells (nx, ny, nz) | spacing [m] | partial volume [m^3] | err % |")
    print("|---|---|---|---|")
    rows = []
    for resolution in ((16, 16, 16), (32, 16, 8), (64, 32, 16)):
        grid = to_voxel_grid(EllipsoidField(semi_axes=ELLIPSOID_SEMI_AXES), resolution)
        volume = grid.volume()
        error = volume / ELLIPSOID_VOLUME - 1.0
        rows.append({"resolution": list(resolution), "spacing": list(grid.spacing), "volume": volume, "error": error})
        print(
            f"| {resolution} | ({grid.spacing[0]:.5f}, {grid.spacing[1]:.5f}, {grid.spacing[2]:.5f}) | "
            f"{volume:.9f} | {100 * error:+.4f} |"
        )
    # Equal cell counts, but only the per-axis grid keeps the spacing near-isotropic.
    isotropic = to_voxel_grid(EllipsoidField(semi_axes=ELLIPSOID_SEMI_AXES), (32, 32, 32))
    anisotropic = to_voxel_grid(EllipsoidField(semi_axes=ELLIPSOID_SEMI_AXES), (64, 32, 16))
    print(
        f"\nSame-body comparison at equal total cost (±3 %): isotropic 32^3 = "
        f"{isotropic.volume() / ELLIPSOID_VOLUME - 1.0:+.4%} on {isotropic.cells:,} cells, "
        f"per-axis (64, 32, 16) = {anisotropic.volume() / ELLIPSOID_VOLUME - 1.0:+.4%} on "
        f"{anisotropic.cells:,} cells.\n"
    )
    assert abs(rows[-1]["error"]) < 1e-3, rows[-1]
    assert anisotropic.spacing[0] < isotropic.spacing[0]
    return {"rows": rows, "isotropic_error": isotropic.volume() / ELLIPSOID_VOLUME - 1.0,
            "anisotropic_error": anisotropic.volume() / ELLIPSOID_VOLUME - 1.0}


def suboff_table() -> dict[str, object]:
    field = SuboffHullField()
    print("## SUBOFF bare hull (AFF-1): the published displacement volume\n")
    print(f"published values: volume {SUBOFF_VOLUME:.6f} m^3 (Groves et al. 1989), "
          f"wetted area {SUBOFF_WETTED_AREA:.6f} m^2 (field quadrature)\n")
    print("| cells (nx, ny, nz) | cell size [m] | partial volume [m^3] | err % | hard-count err % | occupancy | time [s] |")
    print("|---|---|---|---|---|---|---|")
    rows = []
    for resolution in ((109, 13, 13), (218, 26, 26), (436, 52, 52)):
        start = time.perf_counter()
        grid = to_voxel_grid(field, resolution)
        elapsed = time.perf_counter() - start
        partial_error = grid.volume() / SUBOFF_VOLUME - 1.0
        hard_error = grid.occupied_volume() / SUBOFF_VOLUME - 1.0
        rows.append(
            {
                "resolution": list(resolution),
                "cell_size": grid.characteristic_cell_size,
                "volume": grid.volume(),
                "error": partial_error,
                "error_hard": hard_error,
                "inside_fraction": grid.inside_count() / grid.cells,
                "seconds": elapsed,
            }
        )
        print(
            f"| {resolution} | {grid.characteristic_cell_size:.5f} | {grid.volume():.9f} | "
            f"{100 * partial_error:+.4f} | {100 * hard_error:+.4f} | "
            f"{100 * grid.inside_count() / grid.cells:.2f} % | {elapsed:.2f} |"
        )
    assert abs(rows[0]["error"]) < 5e-3, rows[0]
    assert abs(rows[-1]["error"]) < abs(rows[0]["error"]), rows
    # The hard-count error is first order and much larger; that gap is the discretisation
    # error a boolean occupancy grid would hand the solver.
    assert abs(rows[0]["error_hard"]) > abs(rows[0]["error"]), rows[0]
    return {"rows": rows, "published_volume": SUBOFF_VOLUME, "published_area": SUBOFF_WETTED_AREA}


def handoff_record() -> dict[str, object]:
    print("## Solver handoff record\n")
    grid = to_voxel_grid(SuboffHullField(), (218, 26, 26), dtype="float32")
    record = grid.as_handoff()
    print("```json")
    print(json.dumps(record, indent=2)[:1200])
    print("```\n")
    assert record["schema"] == "tensorcad.geometry.voxel_handoff/1"
    assert record["inside_sign"] == "negative_inside"
    assert record["sampling"] == "cell_centred"
    assert record["shape"] == [218, 26, 26]
    assert record["sdf_dtype"] == "float32"
    assert grid.sdf.nbytes == 218 * 26 * 26 * 4
    assert record["sdf_sha256"] == grid.sdf_digest()
    return {"record": record, "nbytes": grid.sdf.nbytes}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None, help="directory for the JSON artifacts")
    arguments = parser.parse_args()

    results = {
        "sphere": sphere_table(),
        "ellipsoid": ellipsoid_table(),
        "suboff": suboff_table(),
        "handoff": handoff_record(),
    }

    if arguments.output is not None:
        arguments.output.mkdir(parents=True, exist_ok=True)
        (arguments.output / "voxel_export.json").write_text(json.dumps(results, indent=2) + "\n")
        (arguments.output / "suboff_voxel_handoff.json").write_text(
            json.dumps(results["handoff"]["record"], indent=2) + "\n"
        )
        print(f"artifacts written to {arguments.output}")

    print("voxel export: all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
