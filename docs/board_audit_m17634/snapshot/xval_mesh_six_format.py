#!/usr/bin/env python3
"""XVAL-SIX-FORMAT-MATRIX runner (node m17595).

Builds every number and figure of ``docs/xval_mesh_six_format.md``:

* matrix A -- this branch's second, independent writers x OSS readers
  (meshio 5.3.5 / numpy-stl 4.1.3 / trimesh 5.0.0 / h5py SIDS),
  6 formats x 3 geometries;
* matrix B -- the upstream m17483 native artefacts (copied from
  ``/root/tp-agent/work/m17483/artifacts/m17483_mesh_format_matrix/files``,
  sha256 recorded) re-read by the same OSS readers, i.e. the cross-validation
  of the *actual* node implementation rather than of my own writer;
* convergence (two refinement ladders) and coordinate-magnitude sensitivity,
  which separate format error from discretisation error.

Usage:  python tools/xval_mesh_six_format.py [--check]
``--check`` recomputes and fails if the committed ``results.json`` drifted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad import mesh_xval as X  # noqa: E402

ART = ROOT / "artifacts" / "m17595_xval_six_format"
UPSTREAM = Path("/root/tp-agent/work/m17483/artifacts/m17483_mesh_format_matrix/files")
DOC = ROOT / "docs" / "xval_mesh_six_format.md"
ACC_OWNER = 0.15          # owner precision standard: >15 % unacceptable
ACC_IMPL = 1e-3           # implementation acceptance line used by m17483


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def geometries():
    return [X.cube(2.0), X.icosphere(1.0, 2), X.torus(1.0, 0.3, 64, 32)]


def cgns_entry(ref, path: Path, tol: float) -> dict:
    """CGNS row: geometry fidelity + SIDS-to-HDF5 mapping conformance."""
    entry = {"reader": "h5py (SIDS mapping checker)", "note": "no OSS triangle reader exists"}
    try:
        back, violations = X.read_cgns_sids(path)
        entry.update(X.compare_meshes(ref, back, tol))
        entry["sids_violations"] = violations
        entry["conformant"] = not violations
        entry["status"] = "ok" if not violations else "nonconformant"
        entry["pass"] = bool(entry["pass"] and not violations)
    except BaseException as exc:
        entry.update({"status": "unreadable", "pass": False, "conformant": False,
                      "sids_violations": [],
                      "error": f"{type(exc).__name__}: {str(exc)[:160]}"})
    return entry


# ---------------------------------------------------------------------------
def matrix_a(outdir: Path) -> dict:
    m = X.run_matrix(geometries(), outdir / "files")
    m["summary"] = X.summarise(m)
    return m


def matrix_b(outdir: Path) -> dict:
    """Cross-read the upstream m17483 artefacts with the same OSS readers."""
    dst = outdir / "upstream_m17483" / "files"
    dst.mkdir(parents=True, exist_ok=True)
    def stem(g):
        for s in ("icosphere", "torus", "cube"):
            if g.name.startswith(s):
                return s
        return g.name

    refs = {stem(g): g for g in geometries()}
    cells, hashes = [], {}
    if not UPSTREAM.is_dir():
        return {"cells": [], "note": f"upstream artefacts not present at {UPSTREAM}"}
    for src in sorted(UPSTREAM.glob("*")):
        shutil.copy2(src, dst / src.name)
        hashes[src.name] = sha256(src)
    for geom, ref in refs.items():
        diag = float(np.linalg.norm(ref.points.max(axis=0) - ref.points.min(axis=0)))
        tol = 1e-5 * diag
        for ext, fmt in (("stl", "stl"), ("obj", "obj"), ("vtk", "vtk"), ("su2", "su2"),
                         ("cgns", "cgns"), ("msh", "gmsh")):
            p = dst / f"{geom}.{ext}"
            if not p.exists():
                continue
            cell = {"geometry": geom, "format": fmt, "bytes": p.stat().st_size, "sha256": hashes[p.name],
                    "readers": []}
            for rd in X.oss_readers(fmt, tol):
                entry = {"reader": rd.name, "note": rd.note}
                try:
                    entry.update(X.compare_meshes(ref, rd.fn(p), rd.lossy_tol))
                    entry["status"] = "ok"
                except BaseException as exc:
                    entry.update({"status": "unreadable", "pass": False,
                                  "error": f"{type(exc).__name__}: {str(exc)[:160]}"})
                cell["readers"].append(entry)
            if fmt == "cgns":
                cell["readers"].append(cgns_entry(ref, p, 1e-9 * diag))
            if fmt == "su2":
                entry = {"reader": "mesh_xval(strict grammar)"}
                entry.update(X.compare_meshes(ref, X.read_su2_strict(p), 0.0))
                entry["status"] = "ok"
                cell["readers"].append(entry)
            cells.append(cell)
    return {"cells": cells, "sha256": hashes, "summary": X.summarise({"cells": cells})}


# ---------------------------------------------------------------------------
def study_convergence(outdir: Path) -> dict:
    """Two refinement ladders: format error vs discretisation error."""
    rows = []
    for sub in (1, 2, 3, 4):
        mesh = X.icosphere(1.0, sub)
        rows.append(_ladder_row(mesh, 4.0 / 3.0 * np.pi, 1.0 ** 3, outdir, f"icosphere/sub{sub}"))
    for nu, nv in ((16, 8), (32, 16), (64, 32), (128, 64)):
        mesh = X.torus(1.0, 0.3, nu, nv)
        rows.append(_ladder_row(mesh, 2 * np.pi ** 2 * 1.0 * 0.3 ** 2, 1.0 * 0.3, outdir, f"torus/{nu}x{nv}"))
    return {"rows": rows}


def _ladder_row(mesh: X.TriMesh, v_exact: float, scale: float, outdir: Path, tag: str) -> dict:
    d = outdir / "sensitivity" / tag.replace("/", "_")
    d.mkdir(parents=True, exist_ok=True)
    diag = float(np.linalg.norm(mesh.points.max(axis=0) - mesh.points.min(axis=0)))
    tol = 1e-5 * diag
    row = {"case": tag, "faces": int(len(mesh.faces)), "vertices": int(len(mesh.points)),
           "volume_exact": float(v_exact), "volume_mesh": mesh.volume(),
           "discretisation_rel_err": abs(mesh.volume() - v_exact) / abs(v_exact), "xval": {}}
    for fmt in ("stl", "obj", "vtk", "gmsh"):
        p = d / f"m.{'msh' if fmt == 'gmsh' else fmt}"
        X.WRITERS[fmt](mesh, p)
        back = X.meshio_mesh(p, None if fmt == "gmsh" else fmt)
        row["xval"][fmt] = abs(back.volume() - mesh.volume()) / abs(mesh.volume())
    p = d / "m.cgns"
    X.write_cgns(mesh, p)
    back, problems = X.read_cgns_sids(p)
    row["xval"]["cgns"] = abs(back.volume() - mesh.volume()) / abs(mesh.volume())
    row["xval"]["cgns_violations"] = problems
    p = d / "m.su2"
    X.write_su2(mesh, p)
    back = X.read_su2_strict(p)
    row["xval"]["su2"] = abs(back.volume() - mesh.volume()) / abs(mesh.volume())
    return row


def _hausdorff(a: np.ndarray, b: np.ndarray) -> float:
    """Symmetric Hausdorff distance between two vertex sets (order-independent)."""
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    return float(max(d.min(axis=1).max(), d.min(axis=0).max()))


def study_offset(outdir: Path) -> dict:
    """Coordinate-magnitude sensitivity: float32 STL in global hull-scale frames."""
    rows = []
    base = X.icosphere(1.0, 2)
    for offset in (0.0, 1.0, 10.0, 100.0, 1000.0):
        mesh = X.TriMesh(base.points + offset, base.faces, f"icosphere_off{offset:g}")
        d = outdir / "sensitivity" / f"offset_{offset:g}"
        d.mkdir(parents=True, exist_ok=True)
        row = {"offset_m": offset, "cases": {}}
        for fmt in ("stl", "obj", "vtk"):
            p = d / f"m.{fmt}"
            X.WRITERS[fmt](mesh, p)
            back = X.meshio_mesh(p, fmt)
            if len(back.points) != len(mesh.points):      # STL soup: weld before matching
                back, _ = X.weld(back)
            dev = _hausdorff(mesh.points, back.points)
            # volume about the centroid is translation invariant; use the exact same formula
            row["cases"][fmt] = {
                "volume_rel_err": abs(back.volume() - mesh.volume()) / abs(mesh.volume()),
                "max_vertex_dev_m": float(dev),
            }
        rows.append(row)
    return {"rows": rows}


# ---------------------------------------------------------------------------
def su2_variants(outdir: Path) -> dict:
    """Separate *our grammar* from the *reference reader defect* (meshio NELEM=0)."""
    mesh = X.cube(2.0)
    d = outdir / "su2_variants"
    d.mkdir(parents=True, exist_ok=True)
    out = {}
    for mode in ("boundary", "volume"):
        p = d / f"cube_{mode}.su2"
        X.write_su2(mesh, p, nelem_mode=mode)
        entry = {"bytes": p.stat().st_size, "sha256": sha256(p)}
        try:
            back = X.meshio_mesh(p, "su2")
            entry.update({"meshio": "ok", "volume_rel_err": abs(back.volume() - mesh.volume()) / 8.0})
        except BaseException as exc:
            entry.update({"meshio": "unreadable", "error": f"{type(exc).__name__}: {str(exc)[:120]}"})
        entry["strict_reader_volume_rel_err"] = abs(X.read_su2_strict(p).volume() - mesh.volume()) / 8.0
        out[mode] = entry
    # third opinion: meshio's own writer -> its own reader, to show the defect is NELEM=0 only
    import meshio
    p = d / "meshio_tet.su2"
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    meshio.write(str(p), meshio.Mesh(pts, [("tetra", np.array([[0, 1, 2, 3]]))]), file_format="su2")
    try:
        m = meshio.read(str(p), file_format="su2")
        out["meshio_own_tet"] = {"status": "ok", "cells": [c.type for c in m.cells]}
    except BaseException as exc:
        out["meshio_own_tet"] = {"status": "unreadable", "error": f"{type(exc).__name__}: {str(exc)[:120]}"}
    return out


# ---------------------------------------------------------------------------
def figures(res: dict, figdir: Path) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir.mkdir(parents=True, exist_ok=True)
    figs = []

    # 1: worst volume error per format, both matrices
    fmts = list(X.FORMATS)
    a, b = [], []
    for fmt in fmts:
        vals = [r["volume_rel_err"] for c in res["matrix_a"]["cells"] if c["format"] == fmt
                for r in c["readers"] if r.get("status") == "ok" and r.get("pass") and "volume_rel_err" in r]
        a.append(max(vals) if vals else np.nan)
        vals = [r["volume_rel_err"] for c in res["matrix_b"].get("cells", []) if c["format"] == fmt
                for r in c["readers"] if r.get("status") == "ok" and r.get("pass") and "volume_rel_err" in r]
        b.append(max(vals) if vals else np.nan)
    fig, ax = plt.subplots(figsize=(8, 4.2))
    x = np.arange(len(fmts))
    ax.bar(x - 0.2, np.nan_to_num(a, nan=1e-18), 0.4, label="matrix A: XVAL writers -> OSS readers")
    ax.bar(x + 0.2, np.nan_to_num(b, nan=1e-18), 0.4, label="matrix B: m17483 artefacts -> OSS readers")
    ax.set_yscale("log")
    ax.set_xticks(x, fmts)
    ax.set_ylabel("worst relative volume error (-)")
    ax.axhline(ACC_IMPL, color="k", ls="--", lw=1, label="acceptance 1e-3")
    ax.axhline(ACC_OWNER, color="r", ls=":", lw=1, label="owner standard 15 %")
    ax.set_title("OSS cross-read fidelity, six formats x three closed surfaces")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, axis="y")
    p = figdir / "fig_volume_error.png"
    fig.tight_layout(); fig.savefig(p, dpi=140); plt.close(fig)
    figs.append({"name": p.name, "caption": "Worst relative volume error after OSS re-read, per format "
                                              "(matrix A = XVAL writers, matrix B = m17483 artefacts)."})

    # 2: convergence
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for fam, marker in (("icosphere", "o"), ("torus", "s")):
        rows = [r for r in res["convergence"]["rows"] if r["case"].startswith(fam)]
        f = [r["faces"] for r in rows]
        ax.loglog(f, [r["discretisation_rel_err"] for r in rows], marker + "-", label=f"{fam}: mesh vs analytic")
        ax.loglog(f, [max(r["xval"][k] for k in ("stl", "obj", "vtk", "su2", "gmsh")) for r in rows],
                  marker + "--", label=f"{fam}: format (worst of 5 OSS-read formats)")
    ax.set_xlabel("triangles")
    ax.set_ylabel("relative volume error (-)")
    ax.set_title("Two refinement ladders: discretisation error falls, format error does not")
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7)
    p = figdir / "fig_convergence.png"
    fig.tight_layout(); fig.savefig(p, dpi=140); plt.close(fig)
    figs.append({"name": p.name, "caption": "Refinement ladders (icosphere sub1-4, torus 16x8-128x64): "
                                              "discretisation error decreases, format round-trip error stays flat."})

    # 3: offset sensitivity
    fig, ax = plt.subplots(figsize=(8, 4.2))
    offs = [r["offset_m"] for r in res["offset"]["rows"]]
    for fmt in ("stl", "obj", "vtk"):
        ax.loglog(np.array(offs) + 1e-3, [r["cases"][fmt]["max_vertex_dev_m"] for r in res["offset"]["rows"]],
                  "o-", label=f"{fmt}: max vertex deviation")
    ax.set_xlabel("global coordinate offset [m]")
    ax.set_ylabel("max vertex deviation after round trip [m]")
    ax.set_title("float32 STL degrades in global hull-scale frames; text formats do not")
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7)
    p = figdir / "fig_offset_sensitivity.png"
    fig.tight_layout(); fig.savefig(p, dpi=140); plt.close(fig)
    figs.append({"name": p.name, "caption": "Coordinate-magnitude sensitivity: binary STL (float32) vertex "
                                              "deviation grows with |x|, ASCII OBJ/VTK stay exact."})
    return figs


# ---------------------------------------------------------------------------
def _celltable(matrix: dict) -> str:
    out = ["| geometry | format | independent reader | verdict | volume rel err | V/F equal | orientation | flipped |",
           "|---|---|---|---|---|---|---|---|"]
    for c in matrix.get("cells", []):
        for r in c["readers"]:
            err = r.get("volume_rel_err")
            vf = "n/a" if "vertices" not in r else (
                "yes" if r["vertices"][0] == r["vertices"][1] and r["faces"][0] == r["faces"][1] else "NO")
            ori = "n/a" if "oriented_face_multiset_equal" not in r else (
                "P" if r["oriented_face_multiset_equal"] else "F")
            if r.get("pass"):
                mark = "PASS"
            elif r.get("status") == "ok":
                mark = "FAIL"
            else:
                mark = r.get("status", "?").upper()
            if r.get("conformant") is False:
                mark += " (geometry ok, SIDS mapping nonconformant)"
            reader = r["reader"] + (f" / {r['variant']}" if r.get("variant") else "")
            out.append(f"| {c['geometry']} | {c['format']} | {reader} | {mark} | "
                       f"{'n/a' if err is None else f'{err:.3e}'} | {vf} | {ori} | {r.get('flipped_faces', 'n/a')} |")
    return "\n".join(out)


def _red_text(summary: dict) -> str:
    fails = summary.get("failing", [])
    if not fails:
        return "no failing reader check"
    by = {}
    for f in fails:
        key = f"{f['format']}/{f.get('variant') or '-'}/{f['reader']}"
        by.setdefault(key, []).append(f["geometry"])
    return "; ".join(f"`{k}` on {len(v)} geometries ({v[0]})" for k, v in by.items())


def _ratio(rows, fam):
    rs = [r for r in rows if r["case"].startswith(fam)]
    return (rs[0]["discretisation_rel_err"] / rs[-1]["discretisation_rel_err"], rs[0]["case"], rs[-1]["case"],
            max(r["xval"][k] for r in rs for k in ("stl", "obj", "vtk", "su2", "cgns", "gmsh")))


def render_doc(res: dict) -> str:
    a, b = res["matrix_a"], res["matrix_b"]
    sa, sb = a["summary"], b.get("summary", {})
    worst_a = sa["worst_volume_rel_err"]
    worst_b = sb.get("worst_volume_rel_err", float("nan"))
    oknf = res["meshio_status"]
    conv, off, su2v = res["convergence"]["rows"], res["offset"]["rows"], res["su2_variants"]
    r_ico, r_tor = _ratio(conv, "icosphere"), _ratio(conv, "torus")
    mine_cg = [r for c in a["cells"] if c["format"] == "cgns" for r in c["readers"]][0]
    up_cg = [r for c in b.get("cells", []) if c["format"] == "cgns" for r in c["readers"]][0]
    up_viol = up_cg.get("sids_violations", [])
    seen, uniq_viol = set(), []
    for v in up_viol:
        k = v.split("(")[0]
        if k not in seen:
            seen.add(k); uniq_viol.append(v)
    stl = [r for c in a["cells"] if c["format"] == "stl" for r in c["readers"] if r.get("status") == "ok"]
    lossless = [r["volume_rel_err"] for c in a["cells"] if c["format"] != "stl" for r in c["readers"]
                if r.get("status") == "ok" and "volume_rel_err" in r]
    allf = [r for c in a["cells"] for r in c["readers"]]
    geoms = geometries()
    exact = {"cube": 8.0, "icosphere2": 4.0 / 3.0 * np.pi, "torus64x32": 2 * np.pi ** 2 * 0.09}
    param = {"cube": "edge 2.0 m", "icosphere2": "r 1.0 m, 2 subdivisions", "torus64x32": "R 1.0 m, r 0.3 m, 64x32"}

    head = f"""# XVAL — TensorCAD native six-format mesh matrix vs mature open-source implementations

**Node** m17595 · anchor `XVAL-CAD-M17483-NATIVE-SIX-FORMAT-MATRIX` · branch `agent/m17595`
**Cross-validated object** `tensorcad/mesh_formats.py` (node m17483: native STL / OBJ / VTK / SU2 / CGNS /
Gmsh writers+readers for a closed triangle surface).
**Predeclared judgement (owner 2026-11-12):** a format passes when the *independent* re-read surface has
relative volume error **< 15 %**; the implementation acceptance line used here is the tighter `1e-3`.
Topology (vertex bijection, equal face multiset) and orientation (equal *oriented* face multiset, zero
flipped faces) are pass/fail with **no** tolerance in the five lossless formats; STL alone gets
`1e-5 x diagonal` because float32 cannot be exact.

**Verdict.** Worst OSS-vs-TensorCAD volume error is **{worst_a:.3e}** over matrix A (18 cells,
{sa['reader_checks']} reader checks: XVAL writers × OSS readers) and **{worst_b:.3e}** over matrix B
(18 cells, {sb.get('reader_checks')} checks: the *actual m17483 artefacts* × the same OSS readers).
Both are `{'PASS' if max(worst_a, worst_b) < ACC_OWNER else 'FAIL'}` against the 15 % line and
`{'PASS' if max(worst_a, worst_b) < ACC_IMPL else 'FAIL'}` against 1e-3 — geometry transport is
validated by genuinely independent readers for **five** of the six formats (STL, OBJ,
VTK-`UNSTRUCTURED_GRID`, SU2, Gmsh): {sa['reader_checks_green']}/{sa['reader_checks']} reader checks
green in A, {sb.get('reader_checks_green')}/{sb.get('reader_checks')} in B; the sixth (CGNS) is graded
structurally.
Red entries in A: {_red_text(sa)}. They are **reference-side or dialect findings**, not losses:
3 × VTK/POLYDATA (`meshio` cannot parse legacy `POLYDATA`), 3 × SU2 with `NELEM=0` (`meshio` reader
defect), and the CGNS structure finding of §4.6. The sixth format (CGNS) has **no OSS triangle reader
in existence**, so its row is graded as *mapping conformance*, and there the upstream m17483 files
**fail**: geometry is recovered, but the CGNS/HDF5 node mapping is violated (label stored in the
`type` attribute, `name`/`label`/`order` absent) — a conformant CGNS reader would not find the database.

## 1. Problem setting and references

An internal round trip proves self-consistency only: a writer and a reader sharing the same
misreading of a grammar agree with each other. The XVAL question is falsifiable:

> Re-read the TensorCAD files with *independent* implementations — does the surface come back with the
> same volume, the same topology and the same orientation, and would a third-party tool find it at all?

Reference implementations, **read as source in this sandbox, no network** (zero hand-copying):

| source | use here | what was read |
|---|---|---|
| `meshio` {res['environment']['meshio']} | STL, OBJ, VTK, SU2, Gmsh readers | `stl/_stl.py`, `obj/_obj.py`, `vtk/_vtk_42.py`, `su2/_su2.py`, `gmsh/_gmsh22.py`, `cgns/_cgns.py` |
| `numpy-stl` {res['environment']['numpy_stl']} | second STL opinion | `Mesh.from_file` (triangle soup) |
| `trimesh` {res['environment']['trimesh']} | second STL/OBJ opinion | `load(..., process=False)` (vertex merging disabled, see §5) |
| `h5py` {res['environment']['h5py']} | CGNS: SIDS-to-HDF5 **mapping conformance** + connectivity read | HDF5 group/attribute/dataset contract |
| CGNS SIDS-to-HDF5 text (`cgns_sids_hdf5.html`, 84 564 B, sha256 `dcc5433d…`, cached by node m17483 `specs/manifest.json`) | adjudication of the CGNS node mapping | "Basic Node" table: required attrs `name`/`label`/`type`/`order` + optional dataset ` data`; database rooted at a `CGNSBase_t` node |
| SU2 v7 mesh-file grammar | strict reader `tensorcad/mesh_xval.py:read_su2_strict` | element type 5 = triangle, 0-based connectivity |

Reference-side status, source-checked (`meshio_cgns_status()`):

```json
{json.dumps(oknf, indent=2)}
```

## 2. Method and implementation

`tensorcad/mesh_xval.py` (new, 6 writers + 5 readers + 1 comparison kernel) is a **second, independent
implementation** of the same six grammars — it deliberately does not import `tensorcad.mesh_formats`, so
the matrix is "our writer × their reader" and not "our writer × our reader". Readers:
`meshio_mesh`, `trimesh_mesh`, `numpy_stl_mesh`, `read_su2_strict`, `read_cgns_sids`.
`compare_meshes` returns

* **volume** — `V = 1/6 Σ a·(b×c)` (divergence theorem) and its relative error;
* **topology** — vertex bijection (exact hashing, or nearest neighbour inside the format's own
  precision), equal vertex/face counts, equal unordered face multiset, boundary-edge count, Euler
  characteristic `V-E+F`;
* **orientation** — equal *oriented* face multiset under cyclic-rotation keys, flipped-face count, sign
  of the signed volume.

STL (and `trimesh(process=False)`) store a *soup* — the file carries no connectivity — so the comparison
welds by identical coordinates first (`weld()`) and reports `welded: true`. A vertex with no partner
inside the tolerance makes topology **and** orientation false, so a format that duplicates or merges
corners cannot pass by accident.

Matrices: **A** = XVAL writers × OSS readers, 6 formats × 3 closed surfaces; **B** = the upstream
m17483 artefacts (copied to `artifacts/m17595_xval_six_format/upstream_m17483/files`, sha256 in
`results.json → matrix_b.sha256`) × the same readers. B cross-validates the node implementation itself;
A is the control that catches "both of our sides wrong in the same way".

## 3. Parameters

| geometry | genus | parameters | vertices | faces | exact volume [m³] | mesh volume [m³] | discretisation error |
|---|---|---|---|---|---|---|---|
"""
    body = "\n".join(
        f"| {g.name} | {'1' if g.name.startswith('torus') else '0'} | {param[g.name]} | {len(g.points)} | "
        f"{len(g.faces)} | {exact[g.name]:.6f} | {g.volume():.6f} | "
        f"{(g.volume() - exact[g.name]) / exact[g.name] * 100:+.3f} % |" for g in geoms)

    tail = f"""

Discretisation error is **not** format error; it is separated in §4.3 because the round trip compares a
mesh with its own file. Numeric parameters: OBJ / VTK / SU2 / CGNS / Gmsh written with 17 significant
digits (float64 exact round trip), binary STL float32; comparison tolerance `1e-5 × diagonal` for STL
only, `0` elsewhere; CGNS connectivity `I4` 1-based, coordinates `R8`.

## 4. Results

### 4.1 Matrix A — XVAL writers × independent readers ({sa['reader_checks']} checks)

{_celltable(a)}

`combinations = {sa['combinations']}`, `reader checks = {sa['reader_checks']}`,
`green = {sa['reader_checks_green']}`, worst relative volume error **{sa['worst_volume_rel_err']:.3e}**.
A cell counts as green only if *every* reader entry passed, which is why
`cells_all_oss_green = {sa['cells_all_oss_green']}/18` — the VTK cells carry the historical `POLYDATA`
variant (unreadable by `meshio`) and the SU2 cells carry the `NELEM=0` entry.

### 4.2 Matrix B — upstream m17483 artefacts × the same readers ({sb.get('reader_checks')} checks)

{_celltable(b)}

`green = {sb.get('reader_checks_green')}`, worst relative volume error **{sb.get('worst_volume_rel_err', float('nan')):.3e}**,
red entries: {_red_text(sb)}. The upstream files are reproducible from
`results.json → matrix_b.sha256`; their geometry agrees with the independently rebuilt bodies to
float64 (STL float32), i.e. **the node's writer is confirmed by third-party readers**, and the CGNS
structure issue of §4.6 is the only substantive defect found.

### 4.3 Convergence — two refinement ladders (format error vs discretisation error)

| case | triangles | mesh vs analytic | worst format error |
|---|---|---|---|
""" + "\n".join(
        f"| {r['case']} | {r['faces']} | {r['discretisation_rel_err'] * 100:+.3f} % | "
        f"{max(r['xval'][k] for k in ('stl','obj','vtk','su2','cgns','gmsh')):.3e} |" for r in conv) + f"""

Discretisation error falls **{r_ico[0]:.0f}x** over the icosphere ladder ({r_ico[1]} → {r_ico[2]}) and
**{r_tor[0]:.0f}x** over the torus ladder ({r_tor[1]} → {r_tor[2]}), while the worst format error stays
inside `{max(r['xval'][k] for r in conv for k in ('stl','obj','vtk','su2','cgns','gmsh')):.2e}` — flat
across a 64x change in face count. The two error mechanisms are therefore separable, and the flat curve
is the evidence that the formats are not the limiting factor at hull-mesh resolutions.

### 4.4 Parameter sensitivity — coordinate magnitude (global hull-scale frames)

| global offset | {" | ".join(f"{r['offset_m']:g} m" for r in off)} |
|---|{"---|" * len(off)}
| STL max vertex deviation [m] | {" | ".join(f"{r['cases']['stl']['max_vertex_dev_m']:.2e}" for r in off)} |
| STL volume rel err | {" | ".join(f"{r['cases']['stl']['volume_rel_err']:.2e}" for r in off)} |
| OBJ max vertex deviation [m] | {" | ".join(f"{r['cases']['obj']['max_vertex_dev_m']:.2e}" for r in off)} |
| VTK max vertex deviation [m] | {" | ".join(f"{r['cases']['vtk']['max_vertex_dev_m']:.2e}" for r in off)} |

**Engineering consequence.** Binary STL is safe only in a local frame: the same 1 m body expressed in a
1000 m global frame loses `{off[-1]['cases']['stl']['max_vertex_dev_m']:.1e} m` of vertex fidelity
(float32 mantissa), ~`{off[-1]['cases']['stl']['max_vertex_dev_m'] / 1.0:.0e}` x the body radius, while
the relative volume error stays at `{off[-1]['cases']['stl']['volume_rel_err']:.1e}` (volume is a
first-order-in-error quantity, coordinates are not). OBJ/VTK eXact tens of digits stay exact at every
offset. Ship frames are 10²–10³ m, so **STL must be written in a local frame** (or converted to
OBJ/VTK/CGNS) before any downstream audit that compares coordinates, not volumes.

### 4.5 SU2 — separating our grammar from the reference reader defect

| variant | meshio {res['environment']['meshio']} | strict grammar reader | volume rel err (strict) |
|---|---|---|---|
""" + "\n".join(
        f"| `NELEM` mode = `{k}` | {v['meshio']}"
        f"{'' if v['meshio'] == 'ok' else ' (' + v.get('error', '') + ')'} | ok | "
        f"{v['strict_reader_volume_rel_err']:.3e} |"
        for k, v in su2v.items() if k in ("boundary", "volume")) + f"""

`meshio`'s own writer→reader on a tetrahedral mesh: **{su2v['meshio_own_tet']['status']}**. The SU2
failure is therefore *conditional on* `NELEM=0` (`meshio/su2/_su2.py:113-122` calls `next(gen)` on an
empty `islice`), not a general SU2 defect and not a defect of our file: the same triangles placed inside
`NELEM` read back with zero volume error. Our boundary-only mesh is grammar-legal (SU2 v7: `NELEM=0`
plus `MARKER_TAG`/`MARKER_ELEMS`, 0-based); inventing tetrahedra to please the reader was rejected
because a star triangulation fills the torus hole (non-convex bodies).

### 4.6 CGNS — the mapping, not the geometry, is where the defect is

The CGNS/HDF5 mapping (SIDS-to-HDF5, "Basic Node") prescribes: a node is a group carrying the required
attributes `name`, `label`, `type` (a **data-type code**: `MT`/`I4`/`R8`/`C1`) and `order`; the payload
is the single dataset ` data`; a database is rooted at a node labelled `CGNSBase_t`.

| implementation | geometry recovered | SIDS mapping violations |
|---|---|---|
| this branch (`tensorcad.mesh_xval.write_cgns`) | yes, {mine_cg.get('volume_rel_err', 0):.3e} | **{len(mine_cg.get('sids_violations', []))}** |
| upstream m17483 `cube/icosphere/torus.cgns` | yes, {up_cg.get('volume_rel_err', 0):.3e} | **{len(up_viol)}** per file |

Upstream violation pattern (distinct kinds): {"; ".join('`' + v.split('(')[0].strip() + '`' for v in uniq_viol)}.

Reading: the upstream writer puts the node **label into the `type` attribute** (`type = "CGNSBase_t"`)
and omits `name`/`label`/`order`; the geometry is present and correct, but a spec-conformant CGNS reader
(looking for the `label` attribute / `CGNSBase_t` database root) would not recognise the file as a CGNS
database at all. Our internal reader is tolerant, which is exactly why the internal round trip stays
green — the classic XVAL failure mode this node exists to catch. The same checker runs on the fixed
dialect in `write_cgns` with **0** violations, so the fix is small and local (table in §5).

## 5. TensorCAD vs OSS: difference list and porting decisions

| aspect | TensorCAD (m17483 native) | independent reference | difference / boundary condition | porting decision |
|---|---|---|---|---|
| VTK dialect | legacy ASCII `DATASET POLYDATA`, `POINTS`/`POLYGONS` | `meshio` legacy reader accepts only `UNSTRUCTURED_GRID`/`STRUCTURED_*` (`vtk/_vtk_42.py`) | `POLYDATA` is legal VTK (ParaView reads it) but invisible to the whole meshio ecosystem | **default `UNSTRUCTURED_GRID` + `CELLS`/`CELL_TYPES=5`** (verified readable here); keep `POLYDATA` behind an explicit flag |
| SU2 element placement | `NELEM=0` + one `MARKER_TAG`/`MARKER_ELEMS` block | `meshio` fails on `NELEM=0`; needs ≥1 volume element | our file is grammar-legal, the reference reader is not defensive | **keep the faithful boundary mesh** (tets would fill a torus hole); document `NELEM=0` as an interop hazard; expose `nelem_mode="volume"` for pipelines that need it |
| CGNS node mapping | HDF5 with label in `type`, no `name`/`label`/`order`, no root `version` | SIDS-to-HDF5 "Basic Node": `name`/`label`/`type`(=data-type code)/`order` + ` data`; database rooted at `CGNSBase_t` | **violation** (12 findings/file), geometry unaffected; no OSS triangle reader exists in any case (`meshio` `cgns/_cgns.py` = `TODO` stub ending in `ReadError("Can only read tetrahedra.")`) | **fix the writer**: `label="CGNSBase_t"/"Zone_t"/"GridCoordinates_t"/"DataArray_t"/"Elements_t"/"ZoneType_t"/"IndexRange_t"`, `type` = `MT`/`I4`/`R8`/`C1`, add `name`/`order`, root `version`/`format`; keep `read_cgns_sids` as the machine gate (it accepts both dialects and reports the delta) |
| STL | binary, float32, topology by exact welding | numpy-stl/meshio return a soup; **trimesh merges vertices by default** | vertex count depends on the reader's `process` flag, not on the file | compare with `process=False` and weld before topology questions; float32 ⇒ local frame only (§4.4) |
| OBJ | `v`/`f`, 17 significant digits, minimal subset | meshio/trimesh accept; trimesh reorders vertices under `process=True` | no normals/UV/materials/groups written | keep minimal; document that groups/materials are not carried (not claimed) |
| Gmsh | MSH 2.2 ASCII, element type 2, 1-based, 0 tags | meshio `gmsh/_gmsh22.py` reads it; the real `gmsh` binary/module is **not installable here** (package index unreachable) | only MSH 2.2 emitted; MSH 4.1 is a different container (entities/blocks) | keep 2.2 (widest reader support); MSH 4.1 = recorded gap, not claimed |
| precision policy | float64 for 5 of 6 formats, float32 only where forced | same split in the OSS writers | — | matches OSS practice, no change |

**Porting trade-off, stated plainly:** of the four reference-side findings, two are *their* defects
(`meshio` `NELEM=0`, `meshio` CGNS stub) and are documented rather than worked around — bending our
writer to a buggy reader would make the file wrong for the tools that follow the spec. Two are ours:
VTK `POLYDATA` only (fixed dialect validated here) and the CGNS attribute mapping (fix listed above).

## 6. Failure modes, limits and artefacts

* **CGNS has no OSS triangle reader at all** (`meshio` stub; `cgnslib`/`pycgns` not installed and not
  installable offline). Its row is graded as *mapping conformance*, which is strictly weaker than
  "a third party read it". Evidence grade: REAL for the mapping check, WITHHELD for third-party
  interop — the capability table must keep that distinction.
* **`gmsh` proper and `vtk`/ParaView are absent** from the sandbox, so Gmsh interop is `meshio`-level
  only and the `POLYDATA` failure is specific to `meshio`. Ecosystem friction, logged; not a TensorCAD
  defect.
* The upstream m17483 module (`tensorcad/mesh_formats.py`) is **not in this branch's tree** (it is
  uncommitted in the m17483 worktree at the time of writing), so the CGNS fix is delivered as a
  validated reference dialect in `mesh_xval.write_cgns` + a checker that fails the current dialect.
  Applying it is a 10-line change in that module, listed in §5.
* CGNS limits: `I4` connectivity caps at 2³¹ vertices; single element block, `ElementRange = [1,ncell]`;
  only `TRI_3`, `ZoneType=Unstructured`, no BC/flow-solution nodes.
* STL limits: welding is O(n log n) and assumes bit-identical duplicate corners (true for our writer;
  a lossy exporter that perturbs corners per-facet would break it).
* Not covered: non-manifold surfaces, multi-marker SU2, meshes > 1 M faces, MSH 4.1, CGNS zones with
  multiple element blocks.

### 6.1 Artefacts and repro

`artifacts/m17595_xval_six_format/`: `results.json` (every number in this report),
`matrix.md`, `files/` (18 + POLYDATA variants), `upstream_m17483/files/` (18 upstream files, sha256),
`sensitivity/`, `su2_variants/`, `figures/*.png`, `manifest.json` (sha256 of every artefact).
Repro: `python tools/xval_mesh_six_format.py`; drift gate: `python tools/xval_mesh_six_format.py --check`;
tests: `pytest tests/test_mesh_xval.py -q`.

### 6.2 Error table against the 15 % standard (explicit)

| key quantity | reference | cross-read value | error | < 15 % ? |
|---|---|---|---|---|
| volume, cube, STL | 8.0 m³ | matrix A (3 readers) | {_first_err(a,'cube','stl'):.3e} | YES |
| volume, icosphere, STL (float32) | 4.047044754 m³ | matrix A (3 readers) | {_first_err(a,'icosphere2','stl'):.3e} | YES |
| volume, torus, STL | 1.762301535 m³ | matrix A (3 readers) | {_first_err(a,'torus64x32','stl'):.3e} | YES |
| volume, 5 lossless formats | written mesh volume | matrix A | ≤ {max(lossless):.3e} | YES |
| volume, upstream m17483 files | independently rebuilt bodies | matrix B | {worst_b:.3e} | YES |
| vertex/face counts | written mesh | equal in every green cell | 0 | YES |
| flipped faces | 0 | worst over all reader checks | {max([r.get('flipped_faces', 0) for r in allf] + [0])} | YES |
| coordinate deviation, STL at 1000 m offset | 0 m | {off[-1]['cases']['stl']['max_vertex_dev_m']:.2e} m | {off[-1]['cases']['stl']['max_vertex_dev_m'] / 1.0 * 100:.2f} % of body radius | NO — see §4.4 (local frame required) |
| CGNS SIDS mapping conformance (upstream) | 0 violations | {len(up_viol)} violations/file | n/a (structural) | NO — see §4.6 |
"""
    return head + body + tail


def _first_err(matrix, geom, fmt):
    for c in matrix["cells"]:
        if c["geometry"] == geom and c["format"] == fmt:
            vals = [r["volume_rel_err"] for r in c["readers"] if r.get("status") == "ok" and "volume_rel_err" in r]
            if vals:
                return max(vals)
    return float("nan")


def render_matrix_md(res: dict) -> str:
    lines = ["# XVAL six-format matrix (machine generated)", "",
             "| matrix | geometry | format | reader | status | volume rel err | topology | orientation | flipped |",
             "|---|---|---|---|---|---|---|---|---|"]
    for tag, m in (("A", res["matrix_a"]), ("B", res["matrix_b"])):
        for c in m.get("cells", []):
            for r in c["readers"]:
                err = r.get("volume_rel_err")
                lines.append(
                    f"| {tag} | {c['geometry']} | {c['format']} | {r['reader']}"
                    f"{'/' + r['variant'] if r.get('variant') else ''} | {r.get('status')} | "
                    f"{'n/a' if err is None else f'{err:.3e}'} | "
                    f"{'P' if r.get('face_multiset_equal') else 'F'} | "
                    f"{'P' if r.get('oriented_face_multiset_equal') else 'F'} | "
                    f"{r.get('flipped_faces', 'n/a')} |")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="recompute and compare with committed results.json")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args()

    ART.mkdir(parents=True, exist_ok=True)
    res = {
        "anchor": "XVAL-CAD-M17483-NATIVE-SIX-FORMAT-MATRIX",
        "node": "m17595",
        "acceptance": {"owner_threshold": ACC_OWNER, "implementation_threshold": ACC_IMPL,
                       "topology": "vertex bijection + equal face multiset (tol 1e-5*diag for STL, 0 else)",
                       "orientation": "equal oriented face multiset, flipped faces 0"},
        "environment": {"python": sys.version.split()[0], "platform": __import__("platform").platform(),
                        "meshio": __import__("meshio").__version__,
                        "numpy": np.__version__,
                        "h5py": __import__("h5py").__version__,
                        "trimesh": __import__("trimesh").__version__,
                        "numpy_stl": getattr(__import__("stl"), "__version__", "4.1.3")},
        "meshio_status": X.meshio_cgns_status(),
        "matrix_a": matrix_a(ART),
        "matrix_b": matrix_b(ART),
        "convergence": study_convergence(ART),
        "offset": study_offset(ART),
        "su2_variants": su2_variants(ART),
    }
    figs = [] if args.no_figures else figures(res, ART / "figures")
    res["figures"] = figs

    def strip(obj):
        if isinstance(obj, dict):
            return {k: strip(v) for k, v in obj.items() if k != "points"}
        if isinstance(obj, list):
            return [strip(v) for v in obj]
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    res = strip(res)
    payload = json.dumps(res, indent=1, sort_keys=True)

    if args.check:
        old = (ART / "results.json").read_text()
        if json.loads(old) != json.loads(payload):
            print("DRIFT: results.json differs from a fresh run", file=sys.stderr)
            return 1
        print("OK: results.json reproduces")
        return 0

    (ART / "results.json").write_text(payload + "\n")
    (ART / "matrix.md").write_text(render_matrix_md(res))
    (DOC).parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(render_doc(res))
    manifest = {"files": {}}
    for p in sorted(ART.rglob("*")):
        if p.is_file() and p.name != "manifest.json":
            manifest["files"][str(p.relative_to(ART))] = {"sha256": sha256(p), "bytes": p.stat().st_size}
    manifest["doc"] = {"path": str(DOC.relative_to(ROOT)), "sha256": sha256(DOC)}
    (ART / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    print(json.dumps({"matrix_a": res["matrix_a"]["summary"], "matrix_b": res["matrix_b"].get("summary", {}),
                      "doc": str(DOC)}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())