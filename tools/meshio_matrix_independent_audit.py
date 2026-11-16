"""Independent audit of the 6-format exchange matrix (anchor GAP-CAD-B2).

The matrix in ``tests/test_meshio_format_matrix.py`` writes every file with the
package and reads it back with the package, so it cannot see a format
misunderstanding the two directions share.  This script is the same 18 cells
(6 formats x cube/sphere/Wigley hull) re-checked from the other side:

**out** (export)
    write with ``tensorcad.meshio_exchange.write_meshio``, then parse the bytes
    on disk with ``tensorcad.format_conformance`` -- spec-level readers that do
    not use meshio and do not use the TensorCAD reader set -- and compare the
    volume / point count / cell count / orientation / closedness the *file*
    declares against the source mesh.

**in** (import)
    write the same geometry with ``meshio`` directly from a plain
    ``meshio.Mesh`` (no TensorCAD writer involved), then read it with
    ``tensorcad.meshio_exchange.read_meshio`` and compare against the source.

Plus, for CGNS, the SIDS conformance of the file (``cgns_issues``): the audit
records how many SIDS violations a meshio-written file has and how many ours
has.

Artefacts (default ``runs/m21004_meshio_matrix18_independent/``):

* ``matrix18.json`` -- every number in the report, one record per cell
* ``files/``        -- the files that were audited (re-checkable without rerunning)
* ``matrix.md``     -- the 18-cell table
* ``figures/``      -- matplotlib panels (volume error, cell-status grid, hull)
* ``manifest.json`` -- sha256 of every artefact

Usage::

    python tools/meshio_matrix_independent_audit.py [--out DIR] [--fig DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import shutil
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

warnings.filterwarnings("ignore", category=RuntimeWarning)

from test_meshio_format_matrix import ANALYTIC_VOLUME, GEOMETRIES  # noqa: E402

from tensorcad import format_conformance as fc  # noqa: E402
from tensorcad.mesh_cgns import write_cgns  # noqa: E402
from tensorcad.meshio_exchange import (  # noqa: E402
    FORMATS_ALL,
    FORMATS_SURFACE,
    FORMATS_VOLUME,
    _EXT,
    TetrahedralMesh,
    read_meshio,
    surface_signed_volume,
    write_meshio,
)

#: The three geometries of the acceptance cell (cube = exact box, sphere =
#: curved body, wigley_hull = the real ship hull with collapsed end stations).
GEOM_NAMES = ("cube", "sphere", "wigley_hull")

#: Volume tolerance per format.  STL stores single-precision coordinates by
#: definition, everything else carries the double values unchanged.
VOLUME_TOL = {"stl": 1e-5}
DEFAULT_TOL = 1e-12


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _rel(path: Path) -> str:
    """Repo-relative path when the file is inside the checkout, else the
    absolute path.  Keeps matrix18.json / manifest.json portable when
    ``--out`` / ``--fig`` point outside the repository (e.g. a scratch dir).
    """
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path.resolve())


def _source_volume(mesh) -> float:
    return mesh.volume() if isinstance(mesh, TetrahedralMesh) else surface_signed_volume(mesh)


def _geom(name: str):
    return next(g for g in GEOMETRIES if g.name == name)


def _relative(new: float, old: float) -> float:
    return abs(new - old) / max(abs(old), 1e-300)


def audit_out(fmt: str, geom_name: str, files_dir: Path) -> dict:
    """Write with TensorCAD, verify by parsing the file bytes independently."""
    geom = _geom(geom_name)
    source = geom.tetra if fmt in FORMATS_VOLUME else geom.mesh
    path = write_meshio(
        files_dir / f"out_{fmt}_{geom_name}{_EXT[fmt]}", source, file_format=fmt
    )
    parsed, meta = fc.parse_file(path, fmt)
    src_vol = _source_volume(source)
    vol = parsed.volume()
    tol = VOLUME_TOL.get(fmt, DEFAULT_TOL)
    record = {
        "direction": "out",
        "format": fmt,
        "geometry": geom_name,
        "topology": "tetra" if isinstance(source, TetrahedralMesh) else "triangle",
        "file": _rel(path),
        "bytes": path.stat().st_size,
        "sha256": _sha256(path),
        "source_volume": src_vol,
        "file_volume": vol,
        "relative_error": _relative(vol, src_vol),
        "volume_tol": tol,
        "source_points": len(source.vertices),
        "file_points": int(len(parsed.points)),
        "source_cells": len(source.tets) if isinstance(source, TetrahedralMesh) else len(source.faces),
        "file_cells": int(len(parsed.cells)),
        "orientation": parsed.orientation(),
        "open_edges": None if parsed.cell_kind == "tetra" else parsed.open_edges(),
        "spec_issues": meta.get("issues", []),
        "encoding": meta.get("encoding"),
    }
    record["volume_ok"] = record["relative_error"] < tol
    record["topology_ok"] = (
        record["file_points"] == record["source_points"]
        and record["file_cells"] == record["source_cells"]
    )
    record["orientation_ok"] = record["orientation"] == "outward"
    record["closed_ok"] = record["open_edges"] in (None, 0)
    record["spec_ok"] = not record["spec_issues"]
    if fmt == "cgns":
        record["cgns_sids_violations"] = fc.cgns_issues(path)
        record["sids_ok"] = not record["cgns_sids_violations"]
    record["green"] = all(
        record[key]
        for key in ("volume_ok", "topology_ok", "orientation_ok", "closed_ok", "spec_ok")
    ) and record.get("sids_ok", True)
    return record


def audit_in(fmt: str, geom_name: str, files_dir: Path) -> dict:
    """Write with meshio alone, then import through the TensorCAD reader."""
    import meshio

    geom = _geom(geom_name)
    if fmt in FORMATS_VOLUME:
        points = np.array(geom.tetra.vertices, dtype=np.float64)
        cells = [("tetra", np.array(geom.tetra.tets, dtype=np.int64))]
        src_vol = geom.tetra.volume()
    else:
        points = np.array(geom.mesh.vertices, dtype=np.float64)
        cells = [("triangle", np.array(geom.mesh.faces, dtype=np.int64))]
        src_vol = surface_signed_volume(geom.mesh)
    path = files_dir / f"in_{fmt}_{geom_name}{_EXT[fmt]}"
    meshio.write(str(path), meshio.Mesh(points, cells), file_format=fmt)
    back = read_meshio(path, file_format=fmt)
    vol = _source_volume(back)
    tol = VOLUME_TOL.get(fmt, DEFAULT_TOL)
    record = {
        "direction": "in",
        "format": fmt,
        "geometry": geom_name,
        "file": _rel(path),
        "sha256": _sha256(path),
        "source_volume": src_vol,
        "file_volume": vol,
        "relative_error": _relative(vol, src_vol),
        "volume_tol": tol,
        "source_points": len(points),
        "file_points": len(back.vertices),
        "source_cells": len(cells[0][1]),
        "file_cells": len(back.tets) if isinstance(back, TetrahedralMesh) else len(back.faces),
    }
    record["volume_ok"] = record["relative_error"] < tol
    record["topology_ok"] = (
        record["file_points"] == record["source_points"]
        and record["file_cells"] == record["source_cells"]
    )
    record["orientation_ok"] = record["file_volume"] > 0
    record["green"] = record["volume_ok"] and record["topology_ok"] and record["orientation_ok"]
    return record


def cgns_defect_evidence(files_dir: Path) -> dict:
    """Document the defect the audit found: meshio's CGNS file is not SIDS."""
    import meshio

    geom = _geom("cube")
    path = write_cgns(files_dir / "cgns_cube_sids.cgns", geom.tetra)
    ours = fc.cgns_issues(path)

    legacy = files_dir / "cgns_cube_meshio_writer.cgns"
    meshio.write(
        str(legacy),
        meshio.Mesh(
            np.array(geom.tetra.vertices, dtype=np.float64),
            [("tetra", np.array(geom.tetra.tets, dtype=np.int64))],
        ),
        file_format="cgns",
    )
    theirs = fc.cgns_issues(legacy)
    # meshio can read our SIDS file (the node names it looks for are kept)...
    meshio_read_ok = True
    try:
        mio = meshio.read(str(path), file_format="cgns")
        meshio_read_ok = mio.cells[0].type == "tetra" and len(mio.cells[0].data) == len(geom.tetra.tets)
    except Exception as exc:  # pragma: no cover - recorded, not raised
        meshio_read_ok = f"{type(exc).__name__}: {exc}"
    return {
        "sids_writer": _rel(path),
        "sids_violations_ours": ours,
        "meshio_writer_file": _rel(legacy),
        "sids_violations_meshio": theirs,
        "meshio_reads_our_sids_file": meshio_read_ok,
    }


def make_figures(rows: list[dict], fig_dir: Path) -> list[dict]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir.mkdir(parents=True, exist_ok=True)
    figures: list[dict] = []

    # 1) volume relative error per cell, log scale
    fig, ax = plt.subplots(figsize=(11, 4.6))
    labels = [f"{r['format']}\n{r['geometry']}" for r in rows]
    values = [max(r["relative_error"], 1e-17) for r in rows]
    colors = ["#2b8a3e" if r["green"] else "#c92a2a" for r in rows]
    ax.bar(range(len(rows)), values, color=colors)
    ax.set_yscale("log")
    ax.axhline(1e-12, color="#495057", ls="--", lw=1, label="tolerance (exact formats)")
    ax.axhline(1e-5, color="#f08c00", ls=":", lw=1, label="tolerance (STL, float32)")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, fontsize=6.5, rotation=90)
    ax.set_ylabel("|V_file - V_source| / V_source")
    ax.set_title("6 formats x 3 geometries: volume fidelity read back from the file bytes")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = fig_dir / "volume_error_matrix.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    figures.append({"name": out.name, "caption": "Relative volume error of each audited cell, read back from the file bytes by the spec parser. Green = inside tolerance and topology/orientation/closedness all pass; red = a failing cell."})

    # 2) status grid, out vs in
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    names = sorted({r["format"] for r in rows})
    geoms = sorted({r["geometry"] for r in rows})
    grid = np.zeros((len(names), len(geoms) * 2))
    for r in rows:
        d = 0 if r["direction"] == "out" else 1
        grid[names.index(r["format"]), geoms.index(r["geometry"]) * 2 + d] = 1.0 if r["green"] else -1.0
    ax.imshow(grid, cmap=matplotlib.colors.ListedColormap(["#c92a2a", "#2b8a3e"]), vmin=-1, vmax=1)
    ax.set_xticks(range(len(geoms) * 2))
    ax.set_xticklabels([f"{g}\n{d}" for g in geoms for d in ("out", "in")], fontsize=8)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names)
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            ax.text(j, i, "pass" if grid[i, j] > 0 else "FAIL", ha="center", va="center", fontsize=7, color="white")
    ax.set_title("18 cells x 2 directions (out = our writer, in = meshio writer)")
    fig.tight_layout()
    out = fig_dir / "cell_status_grid.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    figures.append({"name": out.name, "caption": "Status of all audited cells. 'out' verifies our export by parsing the bytes with an independent spec parser, 'in' verifies our import against files written by meshio alone."})

    # 3) the hull geometry as the file declares it, with the analytic volume
    geom = _geom("wigley_hull")
    parsed, _meta = fc.parse_file(ROOT / _hull_file(rows), "obj")
    tris = parsed.points[parsed.triangles()]
    fig = plt.figure(figsize=(9, 4.4))
    ax = fig.add_subplot(121, projection="3d")
    ax.plot_trisurf(
        tris[:, :, 0].ravel(),
        tris[:, :, 1].ravel(),
        tris[:, :, 2].ravel(),
        triangles=np.arange(len(tris) * 3).reshape(-1, 3),
        color="#4dabf7",
        linewidth=0,
        alpha=0.9,
    )
    ax.set_title("Wigley hull as parsed from the OBJ file")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax2 = fig.add_subplot(122)
    vol_file = parsed.volume()
    analytic = ANALYTIC_VOLUME["wigley_hull"]
    ax2.bar(["file (OBJ)", "analytic"], [vol_file, analytic], color=["#4dabf7", "#adb5bd"])
    ax2.set_ylabel("volume [m^3]")
    ax2.set_title(
        f"hull volume: file vs analytic\nrelative difference {abs(vol_file-analytic)/analytic:.2%}"
    )
    for x, v in enumerate([vol_file, analytic]):
        ax2.text(x, v, f"{v:.6f}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    out = fig_dir / "wigley_hull_from_file.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    figures.append({"name": out.name, "caption": "Wigley hull surface reconstructed from the exported OBJ bytes, and its divergence-theorem volume against the analytic hull volume."})
    return figures


def _hull_file(rows: list[dict]) -> str:
    for r in rows:
        if r["geometry"] == "wigley_hull" and r["format"] == "obj" and r["direction"] == "out":
            return r["file"]
    raise KeyError("no wigley_hull OBJ record")


def markdown_table(rows: list[dict], evidence: dict) -> str:
    lines = [
        "# 6 formats x 3 geometries: independent audit of the exchange matrix",
        "",
        "Anchor: `GAP-CAD-B2-MESH-FORMAT-MATRIX`.",
        "",
        "`out` = written by `tensorcad.meshio_exchange.write_meshio`, verified by parsing the",
        "file bytes with the spec readers in `tensorcad.format_conformance` (no meshio, no",
        "TensorCAD reader). `in` = written by meshio alone from a plain `meshio.Mesh`, read",
        "back by `tensorcad.meshio_exchange.read_meshio`.",
        "",
        "| # | dir | format | geometry | V_file [m3] | V_source [m3] | rel err | pts | cells | orient | closed | spec/SIDS | verdict |",
        "|---|-----|--------|----------|-------------|---------------|---------|-----|-------|--------|--------|-----------|---------|",
    ]
    for i, r in enumerate(rows, 1):
        checks = []
        if r["direction"] == "out":
            checks.append("spec ok" if r["spec_ok"] else f"{len(r['spec_issues'])} spec issue(s)")
            if "sids_ok" in r:
                checks.append("SIDS ok" if r["sids_ok"] else f"{len(r['cgns_sids_violations'])} SIDS violation(s)")
        else:
            checks.append("n/a (meshio writer)")
        lines.append(
            "| {i} | {d} | {f} | {g} | {v:.9f} | {s:.9f} | {e:.2e} | {p} | {c} | {o} | {cl} | {chk} | {verd} |".format(
                i=i,
                d=r["direction"],
                f=r["format"],
                g=r["geometry"],
                v=r["file_volume"],
                s=r["source_volume"],
                e=r["relative_error"],
                p=r["file_points"],
                c=r["file_cells"],
                o="outward" if r["orientation_ok"] else "INWARD",
                cl="-" if r["direction"] == "out" and r.get("open_edges") is None else ("yes" if r.get("closed_ok", True) else "OPEN"),
                chk="; ".join(checks),
                verd="PASS" if r["green"] else "FAIL",
            )
        )
    passing = sum(1 for r in rows if r["green"])
    lines += [
        "",
        f"**{passing}/{len(rows)} audited cells green.**",
        "",
        "## CGNS SIDS conformance",
        "",
        f"* `tensorcad.mesh_cgns.write_cgns`: **{len(evidence['sids_violations_ours'])} SIDS violations**",
        f"* meshio-written CGNS (what the exporter used before this change): "
        f"**{len(evidence['sids_violations_meshio'])} SIDS violations**",
        f"* meshio reads our SIDS file back: {evidence['meshio_reads_our_sids_file']}",
        "",
        "meshio-written file, violations:",
        "",
    ]
    lines += [f"{i}. {p}" for i, p in enumerate(evidence["sids_violations_meshio"], 1)]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        default=str(ROOT / "runs" / "m21004_meshio_matrix18_independent"),
        help="run directory",
    )
    parser.add_argument(
        "--fig",
        default=str(ROOT / "artifacts" / "m21004_meshio_matrix18_independent"),
        help="figure directory",
    )
    args = parser.parse_args(argv)

    run_dir = Path(args.out)
    files_dir = run_dir / "files"
    if files_dir.exists():
        shutil.rmtree(files_dir)
    files_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig) / "figures"

    rows: list[dict] = []
    for fmt in FORMATS_ALL:
        for name in GEOM_NAMES:
            rows.append(audit_out(fmt, name, files_dir))
            rows.append(audit_in(fmt, name, files_dir))

    evidence = cgns_defect_evidence(files_dir)
    figures = make_figures(rows, fig_dir)

    payload = {
        "anchor": "GAP-CAD-B2-MESH-FORMAT-MATRIX",
        "python": platform.python_version(),
        "meshio": __import__("meshio").__version__,
        "geometries": list(GEOM_NAMES),
        "formats": list(FORMATS_ALL),
        "surface_formats": list(FORMATS_SURFACE),
        "volume_formats": list(FORMATS_VOLUME),
        "analytic_volumes": {g: ANALYTIC_VOLUME[g] for g in GEOM_NAMES},
        "rows": rows,
        "cgns": evidence,
        "figures": [_rel(fig_dir / f["name"]) for f in figures],
    }
    (run_dir / "matrix18.json").write_text(json.dumps(payload, indent=2, sort_keys=False))
    (run_dir / "matrix.md").write_text(markdown_table(rows, evidence))

    manifest = {}
    for path in sorted(list(files_dir.glob("*")) + list(fig_dir.glob("*.png")) + [
        run_dir / "matrix18.json",
        run_dir / "matrix.md",
    ]):
        manifest[_rel(path)] = _sha256(path)
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    green = sum(1 for r in rows if r["green"])
    print(f"audited cells: {len(rows)}  green: {green}")
    for r in rows:
        if not r["green"]:
            print("FAIL", r["direction"], r["format"], r["geometry"], r.get("spec_issues"), r.get("cgns_sids_violations"))
    print("CGNS SIDS violations  ours:", len(evidence["sids_violations_ours"]),
          " meshio:", len(evidence["sids_violations_meshio"]))
    print("artefacts:", run_dir, fig_dir)
    return 0 if green == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())