"""Six-format mesh matrix, both directions: round trip *and* foreign intake (m20943).

Anchor: ``GAP-CAD-B2-MESH-FORMAT-MATRIX`` — "meshio 全格式矩阵进出
(STL/OBJ/VTK/SU2/CGNS/Gmsh)", acceptance "往返保真：体积/拓扑/取向不变
（6 格式×3 几何=18 组合全绿表）".

The acceptance table is one axis of this script.  The other axis is the one the
acceptance sentence hides: *进* (intake).  Every round-trip cell writes its file
with our own writer, so 18/18 green says nothing about a quad-dominant cage that
arrives from a CAD kernel, a planar ``NDIME=2`` SU2 file, or a CGNS zone split
over several ``Elements_t`` blocks.  Those are measured here against the same
invariants (volume / topology / orientation) instead of being assumed.

Outputs (``artifacts/m20943_foreign_intake/``):

* ``intake_matrix.json``  every number in the report, programmatically emitted
* ``matrix.md``           markdown tables for the two axes
* ``figures/*.png``       round-trip error chart + intake acceptance matrix
* ``files/``              the foreign fixtures, so a third party can re-run the read

Usage::

    python3 tools/mesh_intake_matrix.py [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

from tensorcad import mesh_formats as mf  # noqa: E402

EXT = {"gmsh": "msh"}
ROUNDTRIP_BODIES = {
    "cube": lambda: mf.cube(),
    "icosphere2": lambda: mf.icosphere(1.0, 2),
    "torus64x32": lambda: mf.torus(),
}
OWNER_THRESHOLD = 0.15
IMPL_THRESHOLD = 1e-3

#: Unit cube as a 6-quad cage (outward winding) — the canonical foreign body.
QUAD_POINTS = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ]
)
QUAD_CAGE = [
    [0, 3, 2, 1],
    [4, 5, 6, 7],
    [0, 1, 5, 4],
    [1, 2, 6, 5],
    [2, 3, 7, 6],
    [3, 0, 4, 7],
]
TRI_CAGE = [face for quad in QUAD_CAGE for face in mf._fan(quad)]


# ---------------------------------------------------------------------------
# foreign fixture writers (deliberately not mf.write: the files must be foreign)
# ---------------------------------------------------------------------------
def write_tri_cage(path: Path, fmt: str) -> None:
    mesh = mf.TriMesh(QUAD_POINTS, np.array(TRI_CAGE, dtype=np.int64), "tri_cage")
    mf.write(mesh, fmt, path)


def write_quad_cage(path: Path, fmt: str) -> None:
    if fmt == "obj":
        path.write_text(
            "".join(f"v {x} {y} {z}\n" for x, y, z in QUAD_POINTS)
            + "".join("f " + " ".join(str(i + 1) for i in q) + "\n" for q in QUAD_CAGE),
            encoding="utf-8",
        )
    elif fmt == "vtk":
        path.write_text(
            "# vtk DataFile Version 3.0\nquad cage\nASCII\nDATASET UNSTRUCTURED_GRID\n"
            + f"POINTS {len(QUAD_POINTS)} float\n"
            + "".join(f"{x} {y} {z}\n" for x, y, z in QUAD_POINTS)
            + f"CELLS {len(QUAD_CAGE)} {5 * len(QUAD_CAGE)}\n"
            + "".join("4 " + " ".join(map(str, q)) + "\n" for q in QUAD_CAGE)
            + f"CELL_TYPES {len(QUAD_CAGE)}\n"
            + "9\n" * len(QUAD_CAGE),
            encoding="utf-8",
        )
    elif fmt == "gmsh":
        path.write_text(
            "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
            + f"$Nodes\n{len(QUAD_POINTS)}\n"
            + "".join(f"{i + 1} {x} {y} {z}\n" for i, (x, y, z) in enumerate(QUAD_POINTS))
            + "$EndNodes\n"
            + f"$Elements\n{len(QUAD_CAGE)}\n"
            + "".join(
                f"{i + 1} 3 0 " + " ".join(str(t + 1) for t in q) + "\n"
                for i, q in enumerate(QUAD_CAGE)
            )
            + "$EndElements\n",
            encoding="utf-8",
        )
    elif fmt == "su2":
        path.write_text(
            "NDIME= 3\n"
            + f"NPOIN= {len(QUAD_POINTS)}\n"
            + "".join(f"{x} {y} {z}\n" for x, y, z in QUAD_POINTS)
            + f"NELEM= {len(QUAD_CAGE)}\n"
            + "".join("9 " + " ".join(map(str, q)) + " 0\n" for q in QUAD_CAGE)
            + "NMARK= 0\n",
            encoding="utf-8",
        )
    elif fmt == "cgns":
        import h5py

        with h5py.File(path, "w") as f:
            mf._sids_attrs(f, {"version": "3.4", "format": "0.9"})
            mf._sids_node(
                f, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4", 1,
                np.array([4.3], dtype=np.float32),
            )
            base = mf._sids_node(f, "Base", "CGNSBase_t", "I4", 1, np.array([3, 3], dtype=np.int32))
            zone = mf._sids_node(
                base, "Zone", "Zone_t", "I4", 1,
                np.array([len(QUAD_POINTS), len(QUAD_CAGE), 0], dtype=np.int32),
            )
            zt = mf._sids_node(zone, "ZoneType", "ZoneType_t", "C1", 1)
            zt.create_dataset(" data", data=np.bytes_("Unstructured")).attrs[" type"] = np.bytes_("C1")
            gc = mf._sids_node(zone, "GridCoordinates", "GridCoordinates_t", "MT", 1)
            for k, axis in enumerate("XYZ"):
                mf._sids_node(
                    gc, f"Coordinate{axis}", "DataArray_t", "R8", k + 1,
                    np.ascontiguousarray(QUAD_POINTS[:, k]),
                )
            el = mf._sids_node(
                zone, "Elements", "Elements_t", "I4", 2, np.array([mf.CGNS_QUAD_4, 0], dtype=np.int32)
            )
            mf._sids_node(
                el, "ElementRange", "IndexRange_t", "I4", 1,
                np.array([1, len(QUAD_CAGE)], dtype=np.int32),
            )
            mf._sids_node(
                el, "ElementConnectivity", "DataArray_t", "I4", 2,
                (np.array(QUAD_CAGE) + 1).astype(np.int32).ravel(),
            )
    else:  # pragma: no cover - STL has no quad primitive by definition
        raise ValueError(f"no quad grammar in {fmt}")


def write_planar_su2(path: Path, _fmt: str | None = None, dim: int = 2) -> None:
    path.write_text(
        f"NDIME= {dim}\nNPOIN= 4\n0.0 0.0\n1.0 0.0\n1.0 1.0\n0.0 1.0\n"
        "NELEM= 2\n5 0 1 2 0\n5 0 2 3 0\nNMARK= 0\n",
        encoding="utf-8",
    )


def write_marker_only_su2(path: Path, _fmt: str | None = None) -> None:
    path.write_text(
        "NDIME= 2\nNPOIN= 4\n0.0 0.0\n1.0 0.0\n1.0 1.0\n0.0 1.0\n"
        "NELEM= 0\nNMARK= 1\nMARKER_TAG= wall\nMARKER_ELEMS= 1\n9 0 1 2 3\n",
        encoding="utf-8",
    )


def write_mixed_su2(path: Path, _fmt: str | None = None) -> None:
    """Triangles and a surface line in one file: the line must be acknowledged."""
    path.write_text(
        "NDIME= 3\nNPOIN= 3\n0 0 0\n1 0 0\n0 1 0\n"
        "NELEM= 1\n5 0 1 2 0\nNMARK= 1\nMARKER_TAG= wall\nMARKER_ELEMS= 1\n3 0 1\n",
        encoding="utf-8",
    )


def write_multiblock_cgns(path: Path, _fmt: str | None = None) -> None:
    import h5py

    tri_block = [[0, 1, 2], [0, 2, 3]]
    quad_block = [[4, 5, 6, 7], [0, 1, 5, 4]]
    with h5py.File(path, "w") as f:
        mf._sids_attrs(f, {"version": "3.4", "format": "0.9"})
        mf._sids_node(f, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4", 1, np.array([4.3], np.float32))
        base = mf._sids_node(f, "Base", "CGNSBase_t", "I4", 1, np.array([3, 3], np.int32))
        zone = mf._sids_node(base, "Zone", "Zone_t", "I4", 1, np.array([8, 4, 0], np.int32))
        zt = mf._sids_node(zone, "ZoneType", "ZoneType_t", "C1", 1)
        zt.create_dataset(" data", data=np.bytes_("Unstructured")).attrs[" type"] = np.bytes_("C1")
        gc = mf._sids_node(zone, "GridCoordinates", "GridCoordinates_t", "MT", 1)
        for k, axis in enumerate("XYZ"):
            mf._sids_node(
                gc, f"Coordinate{axis}", "DataArray_t", "R8", k + 1,
                np.ascontiguousarray(QUAD_POINTS[:, k]),
            )
        first = 1
        for i, (etype, conn) in enumerate(((mf.CGNS_TRI_3, tri_block), (mf.CGNS_QUAD_4, quad_block))):
            el = mf._sids_node(
                zone, f"Elements{i + 1}", "Elements_t", "I4", 2, np.array([etype, 0], np.int32)
            )
            mf._sids_node(
                el, "ElementRange", "IndexRange_t", "I4", 1,
                np.array([first, first + len(conn) - 1], np.int32),
            )
            mf._sids_node(
                el, "ElementConnectivity", "DataArray_t", "I4", 2,
                (np.array(conn) + 1).astype(np.int32).ravel(),
            )
            first += len(conn)


# ---------------------------------------------------------------------------
# measurement
# ---------------------------------------------------------------------------
def measure_roundtrip(outdir: Path) -> list[dict]:
    rows = []
    for fmt in mf.MESH_FORMATS:
        for body_name, build in ROUNDTRIP_BODIES.items():
            body = build()
            path = outdir / "files" / f"rt_{body_name}.{EXT.get(fmt, fmt)}"
            node = {"axis": "roundtrip", "format": fmt, "body": body_name}
            try:
                report = mf.roundtrip(body, fmt, path, tol=1e-5 if fmt == "stl" else 0.0)
                node.update(
                    {
                        "status": "ok" if report["pass"] else "fail",
                        "faces": report["faces"],
                        "vertices": report["vertices"],
                        "volume_ref": report["volume_ref"],
                        "volume_got": report["volume_got"],
                        "volume_rel_err": report["volume_rel_err"],
                        "euler": report["euler"],
                        "boundary_edges": report["boundary_edges"],
                        "face_multiset_equal": report["face_multiset_equal"],
                        "oriented_face_multiset_equal": report["oriented_face_multiset_equal"],
                        "volume_sign_equal": report["volume_sign_equal"],
                    }
                )
            except Exception as exc:  # recorded, never swallowed
                node.update({"status": "error", "error": f"{type(exc).__name__}: {exc}"})
            rows.append(node)
    return rows


def measure_intake(outdir: Path) -> list[dict]:
    cases = {
        "tri_cage": (mf.MESH_FORMATS, write_tri_cage),
        "quad_cage": (("obj", "vtk", "gmsh", "su2", "cgns"), write_quad_cage),
        "planar_ndime2": (("su2",), write_planar_su2),
        "marker_only_quads": (("su2",), write_marker_only_su2),
        "mixed_line_elements": (("su2",), write_mixed_su2),
        "multiblock_cgns": (("cgns",), write_multiblock_cgns),
    }
    rows: list[dict] = []
    for case, (formats, writer) in cases.items():
        for fmt in formats:
            path = outdir / "files" / f"in_{case}.{EXT.get(fmt, fmt)}"
            node = {"axis": "intake", "case": case, "format": fmt}
            try:
                writer(path, fmt)
                notes: list[str] = []
                got = mf.read(path, fmt)
                # STL is a facet soup by definition: weld before counting topology
                topology = mf.weld(got)[0] if fmt == "stl" else got
                node.update(
                    {
                        "status": "ok",
                        "faces": int(len(got.faces)),
                        "welded_for_topology": fmt == "stl",
                        "volume": got.volume(),
                        "area": got.area(),
                        "euler": topology.euler(),
                        "boundary_edges": topology.boundary_edges(),
                        "notes_with_acknowledgement": notes,
                    }
                )
                # the same file must also be refused without notes= when it carries
                # non-surface cells, i.e. silence is never the answer
                if case == "mixed_line_elements":
                    try:
                        mf.read(path, fmt)
                        node["acknowledgement_required"] = False
                    except ValueError:
                        node["acknowledgement_required"] = True
            except ValueError as exc:
                node.update({"status": "refused", "error": str(exc)})
            except Exception as exc:  # pragma: no cover - unexpected
                node.update({"status": "error", "error": f"{type(exc).__name__}: {exc}"})
            rows.append(node)
    return rows


def measure_meshio_parity(outdir: Path) -> list[dict]:
    """Intake of the same foreign file through the independent meshio importer."""
    rows = []
    try:
        import meshio

        from tensorcad.meshio_exchange import meshio_to_mesh
    except Exception as exc:  # pragma: no cover
        return [{"axis": "meshio_parity", "status": "unavailable", "error": str(exc)}]
    for fmt in ("obj", "vtk", "gmsh", "su2"):
        path = outdir / "files" / f"in_quad_cage.{EXT.get(fmt, fmt)}"
        node = {"axis": "meshio_parity", "format": fmt, "case": "quad_cage"}
        try:
            foreign = meshio_to_mesh(meshio.read(str(path), file_format=fmt))
            other = mf.TriMesh(
                np.array(foreign.vertices, dtype=float),
                np.array(foreign.faces, dtype=np.int64),
                "meshio_path",
            )
            native = mf.read(path, fmt)
            node.update(
                {
                    "status": "agree"
                    if {tuple(sorted(f)) for f in other.faces.tolist()}
                    == {tuple(sorted(f)) for f in native.faces.tolist()}
                    else "disagree",
                    "faces_native": int(len(native.faces)),
                    "faces_meshio": int(len(other.faces)),
                    "volume_native": native.volume(),
                    "volume_meshio": other.volume(),
                }
            )
        except Exception as exc:
            node.update({"status": "unreadable", "error": f"{type(exc).__name__}: {exc}"})
        rows.append(node)
    return rows


def make_figures(rt: list[dict], intake: list[dict], figdir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir.mkdir(parents=True, exist_ok=True)
    out = []

    # figure 1: round-trip relative volume error, 6 formats x 3 bodies
    formats = [f for f in mf.MESH_FORMATS if any(r["format"] == f and r["axis"] == "roundtrip" for r in rt)]
    bodies = list(ROUNDTRIP_BODIES)
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    width = 0.8 / len(bodies)
    for k, body in enumerate(bodies):
        vals = []
        for fmt in formats:
            row = next(r for r in rt if r["axis"] == "roundtrip" and r["format"] == fmt and r["body"] == body)
            vals.append(max(row.get("volume_rel_err", np.nan), 1e-18))
        ax.bar([i + k * width for i in range(len(formats))], vals, width, label=body)
    ax.set_yscale("log")
    ax.set_xticks([i + 0.4 - width / 2 for i in range(len(formats))])
    ax.set_xticklabels([f.upper() for f in formats])
    ax.axhline(IMPL_THRESHOLD, ls="--", c="k", lw=1, label=f"implementation gate {IMPL_THRESHOLD:g}")
    ax.axhline(OWNER_THRESHOLD, ls=":", c="r", lw=1, label=f"owner 15% limit")
    ax.set_ylabel("relative volume error after round trip")
    ax.set_title("18-cell round-trip matrix: volume fidelity (6 formats x 3 bodies)")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    p = figdir / "roundtrip_volume_error.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    out.append(p)

    # figure 2: intake acceptance matrix
    case_names = []
    for r in intake:
        if r["case"] not in case_names:
            case_names.append(r["case"])
    formats_i = list(mf.MESH_FORMATS)
    grid = np.full((len(formats_i), len(case_names)), -1.0)
    annotations = {}
    for r in intake:
        i = formats_i.index(r["format"])
        j = case_names.index(r["case"])
        grid[i, j] = {"ok": 1.0, "refused": 0.0}.get(r["status"], 0.5)
        annotations[(i, j)] = {
            "ok": f"ok\n{r.get('faces', '-')}f",
            "refused": "refused",
        }.get(r["status"], r.get("status", "?"))
    fig, ax = plt.subplots(figsize=(9.2, 4.4))
    ax.imshow(grid, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(case_names)))
    ax.set_xticklabels([c.replace("_", "\n") for c in case_names], fontsize=8)
    ax.set_yticks(range(len(formats_i)))
    ax.set_yticklabels([f.upper() for f in formats_i])
    for (i, j), text in annotations.items():
        ax.text(j, i, text, ha="center", va="center", fontsize=7)
    ax.set_title("Foreign-intake matrix: quad cages, planar SU2, marker-only, multi-block CGNS")
    fig.tight_layout()
    p = figdir / "intake_matrix.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    out.append(p)

    # figure 3: what the quad cage became (faces / volume) native vs meshio path
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6))
    labels, native_f, meshio_f = [], [], []
    for r in intake:
        if r["case"] == "quad_cage" and r["status"] == "ok":
            labels.append(r["format"].upper())
            native_f.append(r["faces"])
    parity = {r["format"]: r for r in intake if r["axis"] == "meshio_parity"}
    ax = axes[0]
    ax.bar(labels, native_f, color="tab:blue")
    ax.axhline(2 * len(QUAD_CAGE), ls="--", c="k", lw=1, label="6 quads x 2 = 12")
    ax.set_ylabel("triangles after intake")
    ax.set_title("quad cage -> triangles, native reader")
    ax.legend(fontsize=8)
    ax = axes[1]
    for fmt in formats:
        row = next(r for r in rt if r["axis"] == "roundtrip" and r["format"] == fmt and r["body"] == "cube")
        ax.scatter([fmt.upper()], [abs(row["volume_got"] - 1.0)], s=40)
    ax.set_yscale("log")
    ax.set_ylabel("|V - 1.0| m^3 (cube)")
    ax.set_title("cube volume fidelity")
    fig.tight_layout()
    p = figdir / "intake_faces_and_cube_volume.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    out.append(p)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "artifacts" / "m20943_foreign_intake"))
    args = ap.parse_args()
    outdir = Path(args.out)
    (outdir / "files").mkdir(parents=True, exist_ok=True)

    rt = measure_roundtrip(outdir)
    intake = measure_intake(outdir)
    parity = measure_meshio_parity(outdir)
    figs = make_figures(rt, intake, outdir / "figures")

    payload = {
        "anchor": "GAP-CAD-B2-MESH-FORMAT-MATRIX",
        "mission": 20943,
        "roundtrip": rt,
        "intake": intake,
        "meshio_parity": parity,
        "summary": {
            "roundtrip_cells": len(rt),
            "roundtrip_ok": sum(1 for r in rt if r["status"] == "ok"),
            "intake_cells": len(intake),
            "intake_ok": sum(1 for r in intake if r["status"] == "ok"),
            "intake_refused": sum(1 for r in intake if r["status"] == "refused"),
            "max_roundtrip_rel_err": max(
                (r.get("volume_rel_err", 0.0) for r in rt if r["status"] == "ok"), default=None
            ),
        },
        "figures": [str(p.relative_to(ROOT)) for p in figs],
    }
    (outdir / "intake_matrix.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    lines = ["# Six-format matrix — both directions (m20943)", ""]
    lines += [
        "## Axis A: round trip (the acceptance table)",
        "",
        "| format | body | status | faces (ref/got) | volume rel. err | Euler | boundary edges |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in rt:
        lines.append(
            f"| {r['format']} | {r['body']} | {r['status']} | "
            f"{r.get('faces', ['-', '-'])[0]}/{r.get('faces', ['-', '-'])[1]} | "
            f"{r.get('volume_rel_err', float('nan')):.3e} | {r.get('euler')} | {r.get('boundary_edges')} |"
        )
    lines += [
        "",
        "## Axis B: foreign intake (files TensorCAD did not write)",
        "",
        "| case | format | status | faces | volume / area | Euler | boundary edges |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in intake:
        lines.append(
            f"| {r['case']} | {r['format']} | {r['status']} | {r.get('faces', '-')} | "
            f"{r.get('volume', float('nan')):.6f} / {r.get('area', float('nan')):.6f} | "
            f"{r.get('euler', '-')} | {r.get('boundary_edges', '-')} |"
        )
    lines += ["", "## Cross-check against the meshio importer", "", "| format | status | faces (native/meshio) |", "|---|---|---|"]
    for r in parity:
        lines.append(f"| {r.get('format', '-')} | {r['status']} | {r.get('faces_native', '-')}/{r.get('faces_meshio', '-')} |")
    (outdir / "matrix.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps(payload["summary"], indent=2))
    print("wrote", outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())