"""Run the meshio format matrix (6 formats x 3 geometries) and emit evidence.

Produces, under ``runs/m17802_meshio_matrix18/``:

* ``matrix18.json`` -- one record per (format, geometry) combination with the
  measured volumes, relative error, topology/orientation flags and the sha256 of
  the file that was actually written and read back;
* ``files/`` -- the exported meshes themselves (so the artefacts can be
  re-checked by a third party without re-running the code);
* figures under ``artifacts/m17802_meshio_matrix18/``.

The geometry fixtures are imported from ``tests/test_meshio_format_matrix.py``
on purpose: one definition, and the run therefore exercises exactly the meshes
the test suite gates.

Usage::

    python scripts/meshio_matrix18.py
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

import meshio  # noqa: E402
import numpy as np  # noqa: E402
from test_meshio_format_matrix import ANALYTIC_VOLUME, GEOMETRIES  # noqa: E402

from tensorcad import mesh as meshlib  # noqa: E402
from tensorcad.meshio_exchange import (  # noqa: E402
    _EXT,
    FORMATS_ALL,
    FORMATS_VOLUME,
    TetrahedralMesh,
    round_trip_meshio,
    surface_signed_volume,
    tetrahedralise_surface,
    write_meshio,
)

#: Output locations default to the run that first produced the 18-combination
#: matrix; override with ``--out``/``--fig`` (m18268 added the real hull, i.e.
#: 24 combinations, and writes to its own run directory so the two evidence
#: sets stay separable).
RUN_DIR = ROOT / "runs" / "m17802_meshio_matrix18"
FIG_DIR = ROOT / "artifacts" / "m17802_meshio_matrix18"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _volume(mesh) -> float:
    return mesh.volume() if isinstance(mesh, TetrahedralMesh) else surface_signed_volume(mesh)


def _cells(mesh) -> int:
    return len(mesh.tets) if isinstance(mesh, TetrahedralMesh) else len(mesh.faces)


def run_matrix() -> list[dict]:
    files_dir = RUN_DIR / "files"
    if files_dir.exists():
        shutil.rmtree(files_dir)
    files_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    for fmt in FORMATS_ALL:
        for geom in GEOMETRIES:
            original = geom.tetra if fmt in FORMATS_VOLUME else geom.mesh
            path = write_meshio(files_dir / f"{fmt}_{geom.name}{_EXT[fmt]}", original, file_format=fmt)
            rt = round_trip_meshio(original, file_format=fmt)
            vol_orig, vol_rt = _volume(original), _volume(rt)
            rel = abs(vol_rt - vol_orig) / max(abs(vol_orig), 1e-300)
            limit = 1e-3 if fmt == "stl" else 1e-12
            exact = (
                all(a == b for a, b in zip(original.vertices, rt.vertices, strict=True))
                and (original.tets if isinstance(original, TetrahedralMesh) else original.faces)
                == (rt.tets if isinstance(rt, TetrahedralMesh) else rt.faces)
            )
            rows.append(
                {
                    "format": fmt,
                    "geometry": geom.name,
                    "topology": "tetra" if isinstance(original, TetrahedralMesh) else "triangle",
                    "vertices": len(original.vertices),
                    "cells": _cells(original),
                    "volume_original": vol_orig,
                    "volume_round_trip": vol_rt,
                    "relative_error": rel,
                    "volume_limit": limit,
                    "volume_ok": rel < limit,
                    "topology_ok": len(rt.vertices) == len(original.vertices)
                    and _cells(rt) == _cells(original),
                    "orientation_ok": math.copysign(1.0, vol_rt) == math.copysign(1.0, vol_orig),
                    "edge_on_facets": getattr(original, "edge_on_facets", 0),
                    "exact_restore": bool(exact),
                    "file": str(path.relative_to(ROOT)),
                    "bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                }
            )
    return rows


def geometry_gate() -> list[dict]:
    gate = []
    for geom in GEOMETRIES:
        measured = surface_signed_volume(geom.mesh)
        analytic = ANALYTIC_VOLUME[geom.name]
        tet = tetrahedralise_surface(geom.mesh)
        gate.append(
            {
                "geometry": geom.name,
                "vertices": len(geom.mesh.vertices),
                "faces": len(geom.mesh.faces),
                "watertight": meshlib.is_watertight(geom.mesh),
                "edge_manifold": meshlib.is_edge_manifold(geom.mesh),
                "volume_measured": measured,
                "volume_analytic": analytic,
                "relative_error": abs(measured - analytic) / analytic,
                "tetrahedra": len(tet.tets),
                "tet_volume": tet.signed_volume(),
                "tetV_minus_surfaceV": tet.signed_volume() - measured,
                "edge_on_facets": tet.edge_on_facets,
            }
        )
    return gate


def make_figures(rows: list[dict], gate: list[dict]) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    # Figure 1: volume relative error, 18 combinations, log scale
    labels = [f"{r['format']}\n{r['geometry']}" for r in rows]
    errors = [max(r["relative_error"], 1e-18) for r in rows]
    colours = ["#c0392b" if r["format"] == "stl" else "#1f77b4" for r in rows]
    fig, ax = plt.subplots(figsize=(14, 5))
    bars = ax.bar(range(len(rows)), errors, color=colours)
    ax.set_yscale("log")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("relative volume error after round trip")
    ax.set_title(
        f"meshio format matrix: volume preserved over {len(rows)} combinations "
        f"({len(formats_all := list(dict.fromkeys(r['format'] for r in rows)))} formats x "
        f"{len(list(dict.fromkeys(r['geometry'] for r in rows)))} geometries)"
    )
    ax.axhline(1e-12, color="k", ls="--", lw=1, label="lossless limit 1e-12")
    ax.axhline(1e-3, color="#c0392b", ls=":", lw=1, label="STL limit 1e-3")
    for bar, err in zip(bars, errors, strict=True):
        ax.annotate(f"{err:.1e}", (bar.get_x() + bar.get_width() / 2, err),
                    ha="center", va="bottom", fontsize=6, rotation=90)
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    path = FIG_DIR / "volume_relative_error.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    figures.append(path)

    # Figure 2: status matrix (format x geometry) with volume error annotated
    formats = list(dict.fromkeys(r["format"] for r in rows))
    geoms = list(dict.fromkeys(r["geometry"] for r in rows))
    grid = np.zeros((len(formats), len(geoms)))
    annot = [["" for _ in geoms] for _ in formats]
    for r in rows:
        i, j = formats.index(r["format"]), geoms.index(r["geometry"])
        ok = r["volume_ok"] and r["topology_ok"] and r["orientation_ok"]
        grid[i, j] = 1.0 if ok else 0.0
        annot[i][j] = (
            f"{'PASS' if ok else 'FAIL'}\n"
            f"{'tetra' if r['topology'] == 'tetra' else 'tri'}\n"
            f"{r['vertices']}v/{r['cells']}c\n"
            f"err {r['relative_error']:.1e}"
        )
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.imshow(grid, cmap="RdYlGn", vmin=0, vmax=1)
    ax.set_xticks(range(len(geoms)))
    ax.set_xticklabels(geoms)
    ax.set_yticks(range(len(formats)))
    ax.set_yticklabels(formats)
    for i in range(len(formats)):
        for j in range(len(geoms)):
            ax.text(j, i, annot[i][j], ha="center", va="center", fontsize=7)
    ax.set_title("18 combinations: volume / topology / orientation all preserved")
    fig.tight_layout()
    path = FIG_DIR / "status_matrix.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    figures.append(path)

    # Figure 3: geometry gate -- discrete vs analytic volume
    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(len(gate))
    measured = [g["volume_measured"] for g in gate]
    analytic = [g["volume_analytic"] for g in gate]
    ax.bar(x - 0.2, analytic, 0.4, label="analytic", color="#999999")
    ax.bar(x + 0.2, measured, 0.4, label="discrete surface (divergence theorem)", color="#1f77b4")
    ax.set_xticks(x)
    ax.set_xticklabels([g["geometry"] for g in gate])
    ax.set_ylabel("volume [m^3]")
    for i, g in enumerate(gate):
        ax.annotate(f"{g['relative_error']:.2%}", (i, max(analytic[i], measured[i])),
                    ha="center", va="bottom", fontsize=9)
    ax.legend()
    ax.set_title("geometry gate: fixtures are watertight and within 2 % of analytic volume")
    fig.tight_layout()
    path = FIG_DIR / "geometry_gate.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    figures.append(path)

    return figures


def _parse_outputs(argv: list[str]) -> None:
    """``--out DIR`` / ``--fig DIR`` (both relative to the repo root or absolute)."""
    global RUN_DIR, FIG_DIR
    it = iter(argv)
    for arg in it:
        if arg == "--out":
            RUN_DIR = Path(next(it)).resolve()
        elif arg == "--fig":
            FIG_DIR = Path(next(it)).resolve()
        else:  # pragma: no cover - CLI misuse
            raise SystemExit(f"unknown argument {arg!r}; use --out DIR / --fig DIR")


def main() -> int:
    _parse_outputs(sys.argv[1:])
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    rows = run_matrix()
    gate = geometry_gate()
    figures = make_figures(rows, gate)

    all_ok = all(r["volume_ok"] and r["topology_ok"] and r["orientation_ok"] for r in rows)
    payload = {
        "run": RUN_DIR.name,
        "anchor": "GAP-CAD-B2-MESH-FORMAT-MATRIX",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "meshio": meshio.__version__,
        "numpy": np.__version__,
        "combinations": len(rows),
        "all_combinations_pass": all_ok,
        "geometry_gate": gate,
        "matrix": rows,
        "figures": [str(p.relative_to(ROOT)) for p in figures],
        "reproduce": (
            f"python scripts/meshio_matrix18.py --out {RUN_DIR.relative_to(ROOT)} "
            f"--fig {FIG_DIR.relative_to(ROOT)} && "
            "python -m pytest tests/test_meshio_format_matrix.py -q"
        ),
    }
    (RUN_DIR / "matrix18.json").write_text(json.dumps(payload, indent=2) + "\n")

    print(f"combinations: {len(rows)}  all pass: {all_ok}")
    for r in rows:
        print(
            f"  {r['format']:5s} {r['geometry']:8s} {r['topology']:8s} "
            f"V={r['volume_original']:.9f} rt={r['volume_round_trip']:.9f} "
            f"rel={r['relative_error']:.2e} sha256={r['sha256'][:12]}"
        )
    print("geometry gate:")
    for g in gate:
        print(
            f"  {g['geometry']:8s} watertight={g['watertight']} "
            f"V={g['volume_measured']:.6f} vs {g['volume_analytic']:.6f} "
            f"({g['relative_error']:.2%}), tets={g['tetrahedra']}"
        )
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
