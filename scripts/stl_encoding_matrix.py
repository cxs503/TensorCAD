"""STL encoding matrix (m23297): ASCII vs binary leg, measured, plus the
float32 quantisation envelope as a function of coordinate offset.

Writes, under ``runs/m23297_stl_encoding/``:

* ``matrix.json``   -- one record per (geometry, encoding, offset) with size,
  volume, relative error, max coordinate error, the pre-flight screen and the
  independent (non-meshio) reader's volume;
* ``files/``        -- every STL actually produced, so a third party can
  re-check the artefacts without re-running the writer;
* figures under ``artifacts/m23297_stl_encoding/``.

Usage:  python scripts/stl_encoding_matrix.py
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
sys.path.insert(0, str(ROOT / "tests"))

import meshio  # noqa: E402
import numpy as np  # noqa: E402
from test_meshio_format_matrix import GEOMETRIES  # noqa: E402

from tensorcad.geometry import Mesh  # noqa: E402
from tensorcad.mesh_foreign_dialects import read_stl_binary, write_stl_binary  # noqa: E402
from tensorcad.meshio_exchange import (  # noqa: E402
    EPS32,
    round_trip_meshio,
    stl_quantisation_bound,
    surface_signed_volume,
    write_meshio,
)

RUN_DIR = ROOT / "runs" / "m23297_stl_encoding"
FIG_DIR = ROOT / "artifacts" / "m23297_stl_encoding"
OFFSETS = (0.0, 1.0e1, 1.0e2, 1.0e3, 1.0e4)
ENCODINGS = (("ascii", False), ("binary", True))


def _translate(mesh: Mesh, offset: float) -> Mesh:
    return Mesh(
        vertices=[(x + offset, y, z) for (x, y, z) in mesh.vertices],
        faces=mesh.faces,
        units="m",
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _independent_volume(path: Path) -> tuple[float, float]:
    """Volume and worst stored-normal vs winding-normal deviation, from the
    hand-written spec parser (no meshio involved)."""
    parsed = read_stl_binary(path)
    pts = [tuple(float(c) for c in v) for v in parsed.vertices]
    faces = [(3 * i, 3 * i + 1, 3 * i + 2) for i in range(len(parsed.cells))]
    volume = surface_signed_volume(Mesh(vertices=pts, faces=faces, units="m"))
    worst = 0.0
    normals = parsed.declared_meta["normals"]
    for k, (a, b, c) in enumerate(faces):
        pa, pb, pc = pts[a], pts[b], pts[c]
        u = np.subtract(pb, pa)
        v = np.subtract(pc, pa)
        n = np.cross(u, v)
        norm = float(np.linalg.norm(n))
        if norm == 0.0:
            continue
        n = n / norm
        stored = np.array(normals[k], dtype=float)
        sn = float(np.linalg.norm(stored))
        if sn == 0.0:
            continue
        worst = max(worst, float(np.linalg.norm(n - stored / sn)))
    return volume, worst


def main() -> int:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    (RUN_DIR / "files").mkdir(exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    for case in GEOMETRIES:
        v_ref = surface_signed_volume(case.mesh)
        for encoding, binary in ENCODINGS:
            for offset in OFFSETS:
                body = _translate(case.mesh, offset) if offset else case.mesh
                screen = stl_quantisation_bound(body, encoding=encoding)
                path = RUN_DIR / "files" / f"{case.name}_{encoding}_x{int(offset)}.stl"
                write_meshio(path, body, file_format="stl", stl_binary=binary)
                back = round_trip_meshio(body, file_format="stl", stl_binary=binary)
                v_back = surface_signed_volume(back)
                rel = abs(v_back - surface_signed_volume(body)) / abs(surface_signed_volume(body))
                a = np.sort(np.array(body.vertices, dtype=float), axis=0)
                b = np.sort(np.array(back.vertices, dtype=float), axis=0)
                max_dev = float(np.abs(a - b).max())
                record = {
                    "geometry": case.name,
                    "encoding": encoding,
                    "offset_m": offset,
                    "bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                    "n_vertices": len(back.vertices),
                    "n_faces": len(back.faces),
                    "volume": v_back,
                    "volume_reference": v_ref,
                    "volume_rel_error": rel,
                    "orientation_preserved": v_back * surface_signed_volume(body) > 0,
                    "max_coordinate_error": max_dev,
                    "screen_max_coordinate_error": screen["max_coordinate_error"],
                    "screen_volume_rel_error": screen["screening_volume_rel_error"],
                    "screen_covers_measurement": bool(
                        max_dev <= screen["max_coordinate_error"]
                        and rel <= max(4.0 * screen["screening_volume_rel_error"], 1e-15)
                    ),
                    "screen_vertex_merge_risk": screen["distinct_vertices_may_merge"],
                }
                if binary:
                    ind_volume, normal_dev = _independent_volume(path)
                    record["independent_reader_volume"] = ind_volume
                    record["independent_reader_rel_error"] = abs(ind_volume - v_ref) / abs(v_ref)
                    record["stored_normal_vs_winding_dev"] = normal_dev
                    # Same body through TensorCAD's own spec writer, which
                    # computes the facet normal from the float32 coordinates it
                    # actually stores.
                    own = RUN_DIR / "files" / f"{case.name}_binary_tensorcad_x{int(offset)}.stl"
                    write_stl_binary(own, body, header="tensorcad binary STL")
                    own_volume, own_dev = _independent_volume(own)
                    record["tensorcad_writer_rel_error"] = abs(own_volume - v_ref) / abs(v_ref)
                    record["tensorcad_writer_normal_dev"] = own_dev
                    record["tensorcad_writer_sha256"] = _sha256(own)
                else:
                    # the spec parser refuses ASCII on purpose; record the refusal
                    record["independent_reader_volume"] = None
                    record["independent_reader_rel_error"] = None
                    record["stored_normal_vs_winding_dev"] = None
                    record["tensorcad_writer_rel_error"] = None
                    record["tensorcad_writer_normal_dev"] = None
                    record["tensorcad_writer_sha256"] = None
                rows.append(record)

    payload = {
        "mission": "m23297",
        "anchor": "GAP-CAD-B2-MESH-FORMAT-MATRIX",
        "eps32": EPS32,
        "meshio_version": meshio.__version__,
        "numpy_version": np.__version__,
        "python": platform.python_version(),
        "platform": platform.machine(),
        "offsets_m": list(OFFSETS),
        "geometries": [c.name for c in GEOMETRIES],
        "rows": rows,
    }
    (RUN_DIR / "matrix.json").write_text(json.dumps(payload, indent=2) + "\n")
    # The measurement record is an artefact (tracked in git); the 40 STL payloads
    # stay local and are shipped as one compressed bundle so a third party can
    # re-read exactly the files whose sha256 is in matrix.json.
    import shutil

    shutil.copyfile(RUN_DIR / "matrix.json", FIG_DIR / "matrix.json")
    shutil.make_archive(str(FIG_DIR / "stl_payloads"), "gztar", root_dir=RUN_DIR, base_dir="files")

    _figures(rows)
    _print_table(rows)
    print(f"\nwrote {RUN_DIR/'matrix.json'} ({len(rows)} records)")
    return 0


def _figures(rows: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def get(geom: str, enc: str, off: float, key: str):
        for r in rows:
            if r["geometry"] == geom and r["encoding"] == enc and r["offset_m"] == off:
                return r[key]
        raise KeyError((geom, enc, off, key))

    # Figure 1: relative volume error vs coordinate offset, both encodings.
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for geom, marker in (("sphere", "o"), ("wigley_hull", "s"), ("capsule", "^")):
        xs = [o for o in OFFSETS if o > 0]
        for enc, style in (("binary", "-"), ("ascii", "--")):
            ys = [max(get(geom, enc, o, "volume_rel_error"), 1e-17) for o in xs]
            axes[0].loglog(xs, ys, style, marker=marker, label=f"{geom} {enc}")
    axes[0].set_xlabel("coordinate offset x [m]")
    axes[0].set_ylabel("relative volume error after round trip")
    axes[0].set_title("STL round trip: ASCII is exact, binary float32 is not")
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=7)

    # Figure 2: max coordinate error vs screen bound (binary).
    for geom in ("sphere", "wigley_hull", "capsule"):
        xs = [o if o else 1.0 for o in OFFSETS]
        axes[1].loglog(
            xs,
            [max(get(geom, "binary", o, "max_coordinate_error"), 1e-17) for o in OFFSETS],
            "o-",
            label=f"{geom} measured",
        )
        axes[1].loglog(
            xs,
            [max(get(geom, "binary", o, "screen_max_coordinate_error"), 1e-17) for o in OFFSETS],
            ":",
            label=f"{geom} screen eps32*max|x|",
        )
    axes[1].set_xlabel("coordinate offset x [m] (0 plotted at 1 m)")
    axes[1].set_ylabel("max |Δcoordinate| [m]")
    axes[1].set_title("measured float32 error vs the pre-flight screen")
    axes[1].grid(True, which="both", alpha=0.3)
    axes[1].legend(fontsize=6)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "stl_encoding_envelope.png", dpi=140)
    plt.close(fig)

    # Figure 3: size cost of the ASCII leg.
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    names = [c for c in dict.fromkeys(r["geometry"] for r in rows)]
    ratios = [
        get(n, "ascii", 0.0, "bytes") / get(n, "binary", 0.0, "bytes") for n in names
    ]
    ax.bar(names, ratios, color="#4c72b0")
    ax.axhline(1.0, color="k", lw=0.8)
    ax.set_ylabel("ASCII bytes / binary bytes")
    ax.set_title("interoperability cost of the default (ASCII) STL leg")
    for i, r in enumerate(ratios):
        ax.text(i, r + 0.05, f"{r:.2f}x", ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stl_ascii_cost.png", dpi=140)
    plt.close(fig)

    # Figure 4: stored facet normal vs winding normal (binary), meshio vs TensorCAD writer.
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    for geom, marker in (("sphere", "o"), ("capsule", "s")):
        xs = [o if o else 1.0 for o in OFFSETS]
        ax.loglog(xs, [max(get(geom, "binary", o, "stored_normal_vs_winding_dev"), 1e-12)
                       for o in OFFSETS], marker + "-", label=f"{geom} meshio writer")
        ax.loglog(xs, [max(get(geom, "binary", o, "tensorcad_writer_normal_dev"), 1e-12)
                       for o in OFFSETS], marker + "--", label=f"{geom} TensorCAD writer")
    ax.set_xlabel("coordinate offset x [m] (0 plotted at 1 m)")
    ax.set_ylabel("max |n_stored - n_winding|")
    ax.set_title("binary STL: facet normal computed before vs after the float32 cast")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "stl_normal_consistency.png", dpi=140)
    plt.close(fig)


def _fmt(value, width: int) -> str:
    return "n/a" if value is None else f"{value:{width}.3e}"


def _print_table(rows: list[dict]) -> None:
    head = (
        f"{'geometry':12s} {'enc':7s} {'offset':>9s} {'bytes':>8s} {'relV':>10s} "
        f"{'maxdx':>10s} {'screen':>10s} {'indep relV':>11s} {'meshio ndev':>11s} {'tcad ndev':>10s}"
    )
    print(head)
    print("-" * len(head))
    for r in rows:
        ind = r["independent_reader_rel_error"]
        nd = r["stored_normal_vs_winding_dev"]
        print(
            f"{r['geometry']:12s} {r['encoding']:7s} {r['offset_m']:9.0e} {r['bytes']:8d} "
            f"{r['volume_rel_error']:10.3e} {r['max_coordinate_error']:10.3e} "
            f"{r['screen_max_coordinate_error']:10.3e} "
            f"{'n/a' if ind is None else f'{ind:11.3e}'} "
            f"{'n/a' if nd is None else f'{nd:11.3e}'} "
            f"{_fmt(r['tensorcad_writer_normal_dev'], 10)}"
        )


if __name__ == "__main__":
    raise SystemExit(main())