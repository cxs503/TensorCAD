"""Regenerate the cost/agreement table of ``docs/mesh_quality_numpy_backend.md``.

Run with ``PYTHONPATH=. python scripts/meshquality_numpy_benchmark.py``. It measures the
scalar and the vectorised implementation on the same meshes, checks that every metric of
every cell agrees (including the ``NaN``/``inf`` guard patterns), and writes
``artifacts/mesh_quality/numpy_backend_benchmark.json`` plus a PNG of the per-cell cost.
"""
from __future__ import annotations

import json
import pathlib
import time

import numpy as np

from tensorcad.meshquality.metrics import hex_metrics, tet_metrics, tri_metrics

ROOT = pathlib.Path(__file__).resolve().parents[1]


def grid(n: int, sx: float = 1.0, sy: float = 1.0, sz: float = 0.3):
    xs, ys, zs = np.linspace(0, sx, n + 1), np.linspace(0, sy, n + 1), np.linspace(0, sz, n + 1)
    points = np.stack(np.meshgrid(xs, ys, zs, indexing="ij"), axis=-1).reshape(-1, 3)
    index = lambda i, j, k: (i * (n + 1) + j) * (n + 1) + k  # noqa: E731
    cells = [
        [
            index(i, j, k),
            index(i + 1, j, k),
            index(i + 1, j + 1, k),
            index(i, j + 1, k),
            index(i, j, k + 1),
            index(i + 1, j, k + 1),
            index(i + 1, j + 1, k + 1),
            index(i, j + 1, k + 1),
        ]
        for i in range(n)
        for j in range(n)
        for k in range(n)
    ]
    return points, cells


def cube_tets(n: int):
    points, hexes = grid(n)
    splits = ((0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6))
    tets = [[cell[i] for i in split] for cell in hexes for split in splits]
    return points, tets


def tri_patch(n: int):
    points = [(i / n, j / n, 0.0) for i in range(n + 1) for j in range(n + 1)]
    width = n + 1
    triangles = []
    for i in range(n):
        for j in range(n):
            a = i * width + j
            b, c = a + 1, a + width
            d = c + 1
            triangles += [[a, c, b], [b, c, d]]
    return points, triangles


def compare(scalar, vector) -> dict[str, float | bool]:
    worst = 0.0
    guards = True
    for name in scalar.metrics:
        a = np.asarray(scalar.metrics[name], dtype=float)
        b = np.asarray(vector.metrics[name], dtype=float)
        guards &= bool(np.array_equal(np.isnan(a), np.isnan(b)))
        guards &= bool(np.array_equal(np.isinf(a), np.isinf(b)))
        finite = np.isfinite(a)
        if finite.any():
            worst = max(worst, float(np.max(np.abs(a[finite] - b[finite]))))
    return {"max_abs_diff": worst, "guards_match": guards}


def bench(kind: str, fn, points, cells, rows: list[dict]) -> None:
    start = time.perf_counter()
    scalar = fn(points, cells, backend="scalar")
    scalar_seconds = time.perf_counter() - start
    start = time.perf_counter()
    vector = fn(points, cells, backend="numpy")
    numpy_seconds = time.perf_counter() - start
    row = {
        "kind": kind,
        "cells": len(cells),
        "scalar_s": scalar_seconds,
        "numpy_s": numpy_seconds,
        "us_cell_scalar": 1e6 * scalar_seconds / max(len(cells), 1),
        "us_cell_numpy": 1e6 * numpy_seconds / max(len(cells), 1),
        "speedup": scalar_seconds / numpy_seconds,
    }
    row.update(compare(scalar, vector))
    rows.append(row)
    print(row)


def main() -> int:
    rows: list[dict] = []
    for n in (8, 20, 32):
        points, cells = grid(n)
        bench("hex", hex_metrics, points, cells, rows)
    for n in (4, 8):
        points, cells = cube_tets(n)
        bench("tet", tet_metrics, points, cells, rows)
    for n in (40, 80):
        points, cells = tri_patch(n)
        bench("tri", tri_metrics, points, cells, rows)

    points, cells = grid(100)  # 10^6 hexahedra: the size a CFD case actually has
    start = time.perf_counter()
    hex_metrics(points, cells, backend="numpy")
    million_seconds = time.perf_counter() - start
    print("numpy, 10^6 hexes:", round(million_seconds, 2), "s")

    out = ROOT / "artifacts" / "mesh_quality"
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "rows": rows,
        "numpy_1e6_hex_seconds": million_seconds,
        "env": {"numpy": np.__version__, "host": "loongarch64", "python": "3.11"},
    }
    (out / "numpy_backend_benchmark.json").write_text(json.dumps(payload, indent=2) + "\n")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover - plotting is optional
        return 0
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    labels = [f"{r['kind']}\n{r['cells']} cells" for r in rows]
    x = np.arange(len(rows))
    axes[0].bar(x - 0.2, [r["us_cell_scalar"] for r in rows], 0.4, label="scalar")
    axes[0].bar(x + 0.2, [r["us_cell_numpy"] for r in rows], 0.4, label="numpy")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("microseconds per cell")
    axes[0].set_xticks(x, labels, fontsize=7)
    axes[0].set_title("cost per cell (log scale)")
    axes[0].legend()
    axes[1].bar(x, [r["speedup"] for r in rows], color="tab:green")
    axes[1].set_ylabel("scalar / numpy")
    axes[1].set_xticks(x, labels, fontsize=7)
    axes[1].set_title("speed-up")
    fig.tight_layout()
    fig.savefig(out / "numpy_backend_speedup.png", dpi=130)
    print("wrote", out / "numpy_backend_speedup.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())