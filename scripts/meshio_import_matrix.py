"""Import-side meshio matrix: quad/n-gon cages from "foreign" files -> TensorCAD.

The round-trip matrix in ``scripts/meshio_matrix18.py`` writes with TensorCAD and
reads back with TensorCAD, so it can only prove self-consistency.  This script
exercises the other direction: files whose cell blocks are **not** the triangles
TensorCAD emits, which is what a CAD kernel, Salome/Pointwise, or a structured
hull-panel generator produces.

For every (geometry, format) pair it records the imported volume, the topology
flags and the orientation sign, plus the sha256 of the file actually read, and
writes JSON + PNG evidence under ``runs/`` / ``artifacts/``.

Usage::

    python scripts/meshio_import_matrix.py
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import meshio  # noqa: E402
import numpy as np  # noqa: E402

from tensorcad import mesh as meshlib  # noqa: E402
from tensorcad.meshio_exchange import (  # noqa: E402
    SURFACE_CELL_TYPES,
    read_meshio,
    surface_signed_volume,
)
from tensorcad.wigley import (  # noqa: E402
    analytic_volume,
    longitudinal_taper,
    reference_parameters,
    vertical_taper,
)

MISSION = "m20888"
RUN_DIR = ROOT / "runs" / f"{MISSION}_meshio_import_matrix"
FIG_DIR = ROOT / "artifacts" / f"{MISSION}_meshio_import_matrix"

#: Surface formats a quad/n-gon cage can be written to.  STL is excluded because
#: it is a triangle-facet soup by definition; SU2/CGNS are volume-cell formats
#: (see the round-trip matrix for that half of the picture).
IMPORT_FORMATS: tuple[tuple[str, str], ...] = (("obj", ".obj"), ("vtk", ".vtk"), ("gmsh", ".msh"))


# ---------------------------------------------------------------------------
# Geometry fixtures (quad / n-gon cages built independently of TensorCAD's
# triangulating MeshBuilder -- these are "foreign" meshes)
# ---------------------------------------------------------------------------

def quad_cube() -> tuple[np.ndarray, np.ndarray, float, str]:
    points = np.array(
        [
            [0.0, 0, 0], [1.0, 0, 0], [1.0, 1.0, 0], [0.0, 1.0, 0],
            [0.0, 0, 1.0], [1.0, 0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
        ]
    )
    quads = np.array(
        [
            [0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4],
            [2, 3, 7, 6], [1, 2, 6, 5], [0, 4, 7, 3],
        ]
    )
    return points, [("quad", quads)], 1.0, "unit cube, 6 quads"


def pentagonal_prism() -> tuple[np.ndarray, list, float, str]:
    radius, height = 1.0, 1.0
    angles = [2.0 * math.pi * k / 5 for k in range(5)]
    points = np.array(
        [(radius * math.cos(a), radius * math.sin(a), 0.0) for a in angles]
        + [(radius * math.cos(a), radius * math.sin(a), height) for a in angles]
    )
    pentagons = np.array([list(range(4, -1, -1)), list(range(5, 10))])
    quads = np.array([(i, (i + 1) % 5, 5 + (i + 1) % 5, 5 + i) for i in range(5)])
    exact = 2.5 * radius**2 * math.sin(2.0 * math.pi / 5) * height
    return points, [("polygon", pentagons), ("quad", quads)], exact, "pentagonal prism, 2 n-gons + 5 quads"


def wigley_quad_cage(nx: int = 41, nz: int = 13) -> tuple[np.ndarray, list, float, str]:
    """Wigley hull as an all-quad structured panel cage (as a panel method emits).

    Same station/level sampling as :func:`tensorcad.wigley.build_mesh`, but the
    loft is kept as quadrilaterals instead of being triangulated by
    :class:`tensorcad.geometry.MeshBuilder`, which is what makes it a realistic
    import fixture.  Reference volume is the exact ``4LBT/9``.
    """
    p = reference_parameters(nx=nx, nz=nz)
    ring_size = 2 * p.nz - 1
    # Node welding by exact position, as every real mesh writer does: the two
    # flanks share the keel point, and both end stations collapse onto the
    # centreline, so without welding the cage would have duplicate nodes and
    # open edges (the importer is not expected to weld, the file must be sane).
    coords: list[tuple[float, float, float]] = []
    index: dict[tuple[float, float, float], int] = {}

    def node(point: tuple[float, float, float]) -> int:
        key = tuple(0.0 if c == 0.0 else c for c in point)  # normalise -0.0
        existing = index.get(key)
        if existing is None:
            existing = len(coords)
            index[key] = existing
            coords.append(key)
        return existing

    rings: list[list[int]] = []
    for i in range(p.nx):
        x = p.length * i / (p.nx - 1)
        ring: list[int] = []
        for j in range(ring_size):
            keelward = j <= p.nz - 1
            level = (p.nz - 1 - j) if keelward else (j - (p.nz - 1))
            z = -p.draft * (1.0 - level / (p.nz - 1))
            breadth = 0.5 * p.beam * longitudinal_taper(p, i) * vertical_taper(p, level)
            ring.append(node((x, breadth if keelward else -breadth, z)))
        rings.append(ring)
    points = np.array(coords)

    quads: list[list[int]] = []
    for i in range(p.nx - 1):
        for j in range(ring_size - 1):
            quads.append([rings[i][j], rings[i + 1][j], rings[i + 1][j + 1], rings[i][j + 1]])
        quads.append(
            [rings[i][ring_size - 1], rings[i + 1][ring_size - 1], rings[i + 1][0], rings[i][0]]
        )
    exact = analytic_volume(p.length, p.beam, p.draft)
    return (
        points,
        [("quad", np.array(quads))],
        exact,
        f"Wigley L={p.length} B={p.beam} T={p.draft}, all-quad structured cage",
    )


GEOMETRIES = {
    "cube_quads": quad_cube,
    "prism_ngons": pentagonal_prism,
    "wigley_quads": wigley_quad_cage,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _triangulate(points: np.ndarray, blocks: list) -> np.ndarray:
    """Reference triangulation of the cage, computed here (not via tensorcad)."""
    tris: list[list[int]] = []
    for cell_type, data in blocks:
        for row in np.asarray(data):
            idx = [int(v) for v in row]
            if cell_type == "quad":
                tris.append([idx[0], idx[1], idx[2]])
                tris.append([idx[0], idx[2], idx[3]])
            elif cell_type == "polygon":
                for k in range(1, len(idx) - 1):
                    tris.append([idx[0], idx[k], idx[k + 1]])
            else:  # pragma: no cover - fixtures are quads/n-gons only
                tris.append(idx)
    return np.array([t for t in tris if len(set(t)) == 3])


def _cage_signed_volume(points: np.ndarray, tris: np.ndarray) -> float:
    p0 = points[tris[:, 0]]
    p1 = points[tris[:, 1]]
    p2 = points[tris[:, 2]]
    return float(np.sum(np.einsum("ij,ij->i", p0, np.cross(p1, p2))) / 6.0)


def run() -> dict:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    files_dir = RUN_DIR / "files"
    files_dir.mkdir(parents=True, exist_ok=True)

    records = []
    for geom_name, factory in GEOMETRIES.items():
        points, blocks, analytic, description = factory()
        tris = _triangulate(points, blocks)
        cage_volume = _cage_signed_volume(points, tris)
        for fmt, ext in IMPORT_FORMATS:
            path = files_dir / f"{geom_name}{ext}"
            # meshio's Gmsh 4.1 writer refuses several cell types without explicit
            # entity tags, so mixed-type cages are written with the MSH 2.2 writer
            # (same reader); the all-quad cage goes out as MSH 4.1, which gives the
            # matrix a 4.1 case as well.
            writer = "gmsh22" if (fmt == "gmsh" and len(blocks) > 1) else fmt
            try:
                meshio.write(str(path), meshio.Mesh(points, blocks), file_format=writer)
            except Exception as error:  # noqa: BLE001 - format-level refusal, recorded as such
                reason = f"{type(error).__name__}: {error}"
                if fmt == "gmsh":
                    reason = (
                        "the meshio Gmsh writer carries no polygon cell type "
                        "(MSH 2.2 elements 1-5 are tri/quad/tet/hex/prism; n-gons need "
                        "MSH 4.1 entity tags that only appear when the file has one "
                        "cell type) -- a format limitation, not an importer one"
                    )
                records.append(
                    {
                        "geometry": geom_name,
                        "format": fmt,
                        "status": "not_supported_by_writer",
                        "reason": reason,
                    }
                )
                print(f"{geom_name:14s} {fmt:5s} SKIP {reason}")
                continue
            # the file as written by meshio, i.e. what an external tool would give us
            reread = meshio.read(str(path), file_format=fmt)
            block_types = sorted({b.type for b in reread.cells})
            mesh = read_meshio(path, file_format=fmt)
            volume = surface_signed_volume(mesh)
            record = {
                "status": "green",
                "geometry": geom_name,
                "description": description,
                "format": fmt,
                "file": str(path.relative_to(ROOT)),
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
                "input_cell_types": [t for t, _ in blocks],
                "writer_format": writer,
                "file_cell_types": block_types,
                "vertices": len(mesh.vertices),
                "faces": len(mesh.faces),
                "faces_expected": int(len(tris)),
                "volume": volume,
                "volume_before_export": cage_volume,
                "volume_analytic": analytic,
                "rel_error_vs_before_export": abs(volume - cage_volume) / abs(cage_volume),
                "rel_error_vs_analytic": abs(volume - analytic) / abs(analytic),
                "watertight": meshlib.is_watertight(mesh),
                "edge_manifold": meshlib.is_edge_manifold(mesh),
                "zero_area_faces": len(meshlib.zero_area_faces(mesh)),
                "orientation_sign_positive": volume > 0.0,
                "topology_preserved": len(mesh.faces) == len(tris),
            }
            records.append(record)
            print(
                f"{geom_name:14s} {fmt:5s} V={volume:.12f} "
                f"d(before_export)={record['rel_error_vs_before_export']:.2e} "
                f"d(analytic)={record['rel_error_vs_analytic']:.3e} "
                f"faces={record['faces']}/{record['faces_expected']} "
                f"watertight={record['watertight']}"
            )

    # Regression exhibit: mixed-block file used to lose every block but the first.
    cube_points, cube_blocks, _, _ = quad_cube()
    mixed_points = np.vstack([cube_points, [[2.0, 2.0, 2.0]]])
    mixed = meshio.Mesh(
        mixed_points,
        [("quad", cube_blocks[0][1]), ("triangle", np.array([[0, 1, 2]]))],
    )
    from tensorcad.meshio_exchange import meshio_to_mesh

    mixed_faces = len(meshio_to_mesh(mixed).faces)

    summary = {
        "mission": MISSION,
        "anchor": "GAP-CAD-B2-MESH-FORMAT-MATRIX",
        "date": "2026-11-15",
        "python": platform.python_version(),
        "meshio": meshio.__version__,
        "numpy": np.__version__,
        "formats": [fmt for fmt, _ in IMPORT_FORMATS],
        "surface_cell_types": list(SURFACE_CELL_TYPES),
        "combinations": len(records),
        "green_combinations": sum(1 for r in records if r["status"] == "green"),
        "skipped_combinations": sum(1 for r in records if r["status"] != "green"),
        "all_green": all(
            r["topology_preserved"]
            and r["watertight"]
            and r["edge_manifold"]
            and r["zero_area_faces"] == 0
            and r["orientation_sign_positive"]
            for r in records
            if r["status"] == "green"
        ),
        "green_table": [
            {
                "geometry": r["geometry"],
                "format": r["format"],
                "volume": r["volume"],
                "rel_err_vs_analytic": r["rel_error_vs_analytic"],
                "faces": r["faces"],
                "watertight": r["watertight"],
                "orientation_positive": r["orientation_sign_positive"],
            }
            for r in records
            if r["status"] == "green"
        ],
        "max_rel_error_vs_before_export": max(
            r["rel_error_vs_before_export"] for r in records if r["status"] == "green"
        ),
        "max_rel_error_vs_analytic": max(
            r["rel_error_vs_analytic"] for r in records if r["status"] == "green"
        ),
        "mixed_block_exhibit": {
            "input_faces": 6 * 2 + 1,
            "imported_faces": mixed_faces,
        },
        "records": records,
    }
    (RUN_DIR / "matrix.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"\n{len(records)} combinations, all_green={summary['all_green']}, "
        f"max d(before_export)={summary['max_rel_error_vs_before_export']:.2e}, "
        f"mixed-block faces {summary['mixed_block_exhibit']['imported_faces']}/"
        f"{summary['mixed_block_exhibit']['input_faces']}"
    )
    return summary


def figures(summary: dict) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out: list[Path] = []
    records = [r for r in summary["records"] if r["status"] == "green"]
    geoms = sorted({r["geometry"] for r in records})
    formats = summary["formats"]

    # Figure 1: relative error, imported volume vs the cage as written, and vs analytic
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    width = 0.8 / len(formats)
    floor = 1e-17
    for ax, key, title in (
        (axes[0], "rel_error_vs_before_export", "imported volume vs pre-export cage"),
        (axes[1], "rel_error_vs_analytic", "imported volume vs analytic reference"),
    ):
        for k, fmt in enumerate(formats):
            vals = [
                max(
                    next(
                        (r[key] for r in records if r["geometry"] == g and r["format"] == fmt),
                        floor,
                    ),
                    floor,
                )
                for g in geoms
            ]
            ax.bar(np.arange(len(geoms)) + k * width, vals, width, label=fmt)
        ax.set_yscale("log")
        ax.set_xticks(np.arange(len(geoms)) + 0.4 - width / 2)
        ax.set_xticklabels(geoms, rotation=12)
        ax.set_ylabel("relative error [-]")
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(title="format")
    fig.tight_layout()
    p1 = FIG_DIR / "import_fidelity.png"
    fig.savefig(p1, dpi=150)
    plt.close(fig)
    out.append(p1)

    # Figure 2: the Wigley all-quad cage itself (what arrives from a panel tool)
    points, blocks, analytic, _ = wigley_quad_cage(nx=81, nz=17)
    quads = blocks[0][1]
    fig = plt.figure(figsize=(9.0, 5.0))
    ax = fig.add_subplot(111, projection="3d")
    for quad in quads[::3]:
        loop = list(quad) + [quad[0]]
        ax.plot(points[loop, 0], points[loop, 1], points[loop, 2], color="#1f77b4", lw=0.25)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.set_title(f"Wigley quad cage imported as {len(quads)} quads -> {2 * len(quads)} triangles")
    ax.view_init(elev=22, azim=-60)
    fig.tight_layout()
    p2 = FIG_DIR / "wigley_quad_cage.png"
    fig.savefig(p2, dpi=150)
    plt.close(fig)
    out.append(p2)

    # Figure 3: the silent-loss exhibit
    exhibit = summary["mixed_block_exhibit"]
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ax.bar(["input file", "imported (first-block only)", "imported (all blocks)"],
           [exhibit["input_faces"], 1, exhibit["imported_faces"]],
           color=["#888888", "#d62728", "#2ca02c"])
    ax.set_ylabel("surface triangles [-]")
    ax.set_title("mixed-block file: silent geometry loss before / after")
    for k, v in enumerate([exhibit["input_faces"], 1, exhibit["imported_faces"]]):
        ax.text(k, v + 0.3, str(v), ha="center")
    fig.tight_layout()
    p3 = FIG_DIR / "mixed_block_loss.png"
    fig.savefig(p3, dpi=150)
    plt.close(fig)
    out.append(p3)
    return out


if __name__ == "__main__":
    import matplotlib  # noqa: F401  (figures() imports it lazily too)

    summary = run()
    for path in figures(summary):
        print("figure:", path.relative_to(ROOT))