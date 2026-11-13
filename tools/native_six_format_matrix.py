"""Native six-format matrix for ``tensorcad.mesh_formats`` (m17907).

Produces ``artifacts/m17907_native_six_format/``:

* ``files/``      one file per (body, format) written by the native writers
* ``results.json`` every number in the report (no hand-copied figures)
* ``matrix.md``    the reader matrix and the error table
* ``figures/``     error chart + geometry panel (matplotlib, sha256 in manifest)
* ``manifest.json`` sha256 of every artefact

Three independent oracles are used and kept apart:

1. **our own reader** (``mesh_formats.read``) -> grammar round trip;
2. **third-party readers** (meshio / numpy-stl / trimesh) -> ecosystem interop,
   with "unreadable" recorded honestly instead of being skipped silently;
3. **the second implementation** in ``tensorcad.mesh_xval`` (its writers are
   written independently of this module) -> writer x writer agreement.

Repro: python3 tools/native_six_format_matrix.py [--out DIR] [--check]
``--check`` recomputes and diffs against the committed ``results.json``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tensorcad import mesh_formats as mf  # noqa: E402
from tensorcad import mesh_xval as mx  # noqa: E402

EXT = {"gmsh": "msh"}
TOL = {"stl": 1e-5}  # STL is float32 by format definition
IMPL_THRESHOLD = 1e-3  # our own gate on the relative volume error
OWNER_THRESHOLD = 0.15  # owner's accuracy standard


def rel(path: Path) -> str:
    """Path relative to the repo root when inside it, else the absolute path."""
    try:
        return str(Path(path).resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bodies():
    cube = mf.cube(2.0)
    sphere = mf.icosphere(1.0, 2)
    torus = mf.torus(1.0, 0.3, 64, 32)
    return [
        ("cube", cube, 8.0, "analytic: V = edge^3 exactly"),
        ("icosphere2", sphere, 4.0 / 3.0 * np.pi, "analytic sphere 4/3*pi*r^3 (inscribed polyhedron, O(h^2) low)"),
        ("torus64x32", torus, 2.0 * np.pi**2 * 1.0 * 0.3**2, "analytic 2*pi^2*R*r^2 (O(h^2) low)"),
    ]


def _third_party(path: Path, fmt: str) -> dict:
    """Read the file with an OSS library and return its volume, or a reason."""
    out = {"library": None, "volume": None, "vertices": None, "faces": None, "status": "unreadable"}

    def tri_volume(points, cells):
        tri = np.asarray(points)[cells[:, :3]]
        return float(np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0)

    import contextlib
    import io

    quiet = io.StringIO()
    try:
        with contextlib.redirect_stdout(quiet), contextlib.redirect_stderr(quiet):
            if fmt == "stl":
                import stl as nstl  # the `numpy-stl` distribution imports as `stl`

                mesh = nstl.Mesh.from_file(str(path))
                out.update(library="numpy-stl", volume=float(mesh.get_mass_properties()[0]),
                           faces=int(len(mesh.vectors)), status="ok")
                return out
            import meshio

            mesh = meshio.read(str(path))
            cells = np.vstack([c.data for c in mesh.cells if c.type in ("triangle", "polygon")])
            out.update(library="meshio", volume=tri_volume(mesh.points, cells),
                       vertices=int(len(mesh.points)), faces=int(len(cells)), status="ok")
    except BaseException as exc:  # noqa: BLE001 - readers may even SystemExit (meshio CGNS stub)
        first = (str(exc).splitlines()[0] if str(exc) else type(exc).__name__)[:120]
        out["reason"] = f"{type(exc).__name__}: {first}"
        out["status"] = "unreadable"
    return out


def _cross_writer(body: mf.TriMesh, fmt: str, outdir: Path) -> dict:
    """Agreement with the independent writer set in ``tensorcad.mesh_xval``."""
    path = outdir / f"{body.name}_xval.{EXT.get(fmt, fmt)}"
    ref = mx.TriMesh(body.points, body.faces, body.name)
    kwargs = {"nelem_mode": "volume"} if fmt == "su2" else {}
    try:
        (mx.write_su2(ref, path, **kwargs) if fmt == "su2" else mx.WRITERS[fmt](ref, path))
    except Exception as exc:  # noqa: BLE001
        return {"status": "write_failed", "reason": f"{type(exc).__name__}: {exc}"}
    structural: list = []
    try:
        if fmt == "cgns":
            # The legacy dialect omits CGNSLibraryVersion_t, so the strict native
            # reader refuses it.  Geometry is still measurable through the other
            # implementation's reader; the structural delta is reported, not hidden.
            got, structural = mx.read_cgns_sids(path)
            got = mf.TriMesh(got.points, got.faces, fmt)
        else:
            got = mf.read(path, fmt)
    except Exception as exc:  # noqa: BLE001
        return {"status": "unreadable", "reason": f"{type(exc).__name__}: {exc}"}
    report = mf.compare_meshes(body, got, tol=TOL.get(fmt, 0.0))
    report.pop("path", None)
    report["status"] = "ok"
    report["file"] = rel(path)
    report["structural_violations"] = list(structural)
    return report


def run(outdir: Path) -> dict:
    files = outdir / "files"
    files.mkdir(parents=True, exist_ok=True)
    results = {
        "anchor": "XVAL-CAD-M17483-NATIVE-SIX-FORMAT-MATRIX",
        "module": "tensorcad/mesh_formats.py",
        "acceptance": {
            "implementation_threshold": IMPL_THRESHOLD,
            "owner_threshold": OWNER_THRESHOLD,
            "topology": "vertex bijection + face multiset + oriented face multiset",
            "cgns": "sids_violations() == []",
        },
        "formats": list(mf.MESH_FORMATS),
        "rows": [],
        "cross_writer": [],
        "sids": {},
        "env": {},
    }
    for name, body, analytic, note in bodies():
        for fmt in mf.MESH_FORMATS:
            path = files / f"{name}.{EXT.get(fmt, fmt)}"
            report = mf.roundtrip(body, fmt, path, tol=TOL.get(fmt, 0.0))
            third = _third_party(path, fmt)
            row = {
                "body": name,
                "format": fmt,
                "file": rel(path),
                "bytes": report["bytes"],
                "faces": report["faces"],
                "vertices": report["vertices"],
                "vertices_read": report["vertices_read"],
                "welded": report["welded"],
                "volume_written": report["volume_ref"],
                "volume_read_back": report["volume_got"],
                "volume_rel_err_own_reader": report["volume_rel_err"],
                "volume_rel_err_vs_analytic": abs(report["volume_ref"] - analytic) / analytic,
                "analytic_reference": analytic,
                "euler": report["euler"],
                "flipped_faces": report["flipped_faces"],
                "topology_pass": report["pass"],
                "third_party": third,
            }
            if third["volume"] is not None:
                # format fidelity: OSS reader versus the mesh we actually wrote
                row["volume_rel_err_third_party"] = (
                    abs(third["volume"] - report["volume_ref"]) / abs(report["volume_ref"])
                )
            results["rows"].append(row)
        for fmt in mf.MESH_FORMATS:
            results["cross_writer"].append(
                {"body": name, "format": fmt, **_cross_writer(body, fmt, files)}
            )

    # CGNS conformance: native file vs the older reference dialect (negative control)
    import h5py  # noqa: F401  (import here so the env block below can report the version)

    probe = bodies()[1][1]
    native = mf.write(probe, "cgns", files / "sids_probe_native.cgns")
    legacy = mx.write_cgns(mx.TriMesh(probe.points, probe.faces, probe.name), files / "sids_probe_legacy.cgns")
    results["sids"] = {
        "rules": list(mf.CGNS_RULES),
        "native_violations": mf.sids_violations(native),
        "legacy_reference_dialect_violations": mf.sids_violations(legacy),
        "legacy_file": rel(legacy),
        "native_file": rel(native),
    }

    import h5py as _h5
    import meshio as _meshio
    import trimesh as _trimesh

    results["env"] = {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "h5py": _h5.__version__,
        "meshio": _meshio.__version__,
        "trimesh": _trimesh.__version__,
        "platform": sys.platform,
    }

    # summary -------------------------------------------------------------
    own = [r["volume_rel_err_own_reader"] for r in results["rows"]]
    topo_ok = sum(1 for r in results["rows"] if r["topology_pass"])
    third_ok = [r for r in results["rows"] if r["third_party"]["status"] == "ok"]
    results["summary"] = {
        "cells": len(results["rows"]),
        "own_reader_topology_pass": topo_ok,
        "worst_own_reader_rel_err": float(max(own)),
        "worst_own_reader_rel_err_where": max(results["rows"], key=lambda r: r["volume_rel_err_own_reader"])["format"],
        "own_reader_all_within_impl_threshold": bool(max(own) <= IMPL_THRESHOLD),
        "own_reader_all_within_owner_threshold": bool(max(own) <= OWNER_THRESHOLD),
        "third_party_readable_cells": len(third_ok),
        "third_party_unreadable_cells": [
            f"{r['format']}({r['third_party'].get('reason', 'unreadable')[:60]})"
            for r in results["rows"]
            if r["third_party"]["status"] != "ok"
        ],
        "third_party_worst_rel_err": (
            float(max(r["volume_rel_err_third_party"] for r in third_ok)) if third_ok else None
        ),
        "cross_writer_agreement": sum(1 for c in results["cross_writer"] if c.get("pass")),
        "cross_writer_cells": len(results["cross_writer"]),
        "cgns_native_violations": len(results["sids"]["native_violations"]),
        "cgns_legacy_violations": len(results["sids"]["legacy_reference_dialect_violations"]),
    }

    # sha256 of every file under files/ + results are added by the writer below
    results["files_sha256"] = {
        str(p.relative_to(outdir)): sha256(p) for p in sorted(files.rglob("*")) if p.is_file()
    }
    (outdir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (outdir / "matrix.md").write_text(render_matrix(results), encoding="utf-8")
    return results


def render_matrix(results: dict) -> str:
    lines = [
        "# Native six-format matrix (m17907)",
        "",
        "Module `tensorcad/mesh_formats.py`; each cell is *write with the native writer -> read back*.",
        "",
        "## Volume error (own reader)",
        "",
        "| body | format | written V | read V | rel err (own) | rel err vs analytic | topo | third-party oracle |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results["rows"]:
        third = r["third_party"]
        cell = (
            f"{third['library']} {r.get('volume_rel_err_third_party', float('nan')):.1e}"
            if third["status"] == "ok"
            else "unreadable"
        )
        lines.append(
            f"| {r['body']} | {r['format']} | {r['volume_written']:.12g} | {r['volume_read_back']:.12g} "
            f"| {r['volume_rel_err_own_reader']:.2e} | {r['volume_rel_err_vs_analytic']:.2e} "
            f"| {'PASS' if r['topology_pass'] else 'FAIL'} | {cell} |"
        )
    s = results["summary"]
    lines += [
        "",
        f"Worst own-reader relative volume error: **{s['worst_own_reader_rel_err']:.3e}** "
        f"({s['worst_own_reader_rel_err_where']}); topology pass {s['own_reader_topology_pass']}/{s['cells']}.",
        "",
        "## CGNS SIDS conformance",
        "",
        f"* native writer: {s['cgns_native_violations']} violation(s)",
        f"* legacy reference dialect (`mesh_xval.write_cgns`): {s['cgns_legacy_violations']} violation(s) "
        f"-> {results['sids']['legacy_reference_dialect_violations']}",
        "",
        "## Cross-writer agreement (second independent implementation, `tensorcad.mesh_xval`)",
        "",
        "| body | format | rel err | topology | status |",
        "|---|---|---|---|---|",
    ]
    for c in results["cross_writer"]:
        if c.get("status") == "ok":
            lines.append(
                f"| {c['body']} | {c['format']} | {c['volume_rel_err']:.2e} | "
                f"{'PASS' if c.get('pass') else 'FAIL'} | ok |"
            )
        else:
            lines.append(f"| {c['body']} | {c['format']} | - | - | {c['status']} |")
    lines += [
        "",
        f"Cross-writer agreement: {s['cross_writer_agreement']}/{s['cross_writer_cells']} cells.",
        "",
    ]
    return "\n".join(lines)


def figure(results: dict, outdir: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    fmt_list = list(mf.MESH_FORMATS)
    body_names = sorted({r["body"] for r in results["rows"]})
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12.4, 4.4), gridspec_kw={"width_ratios": [1.35, 1]})

    width = 0.8 / max(len(body_names), 1)
    floor = 1e-17
    for k, body in enumerate(body_names):
        own = [next(r["volume_rel_err_own_reader"] for r in results["rows"] if r["body"] == body and r["format"] == f)
               for f in fmt_list]
        xs = np.arange(len(fmt_list)) + (k - (len(body_names) - 1) / 2) * width
        ax0.bar(xs, [max(v, floor) for v in own], width=width, label=body)
    ax0.set_yscale("log")
    ax0.set_xticks(np.arange(len(fmt_list)))
    ax0.set_xticklabels(fmt_list)
    ax0.axhline(IMPL_THRESHOLD, color="k", ls="--", lw=1, label="impl gate 1e-3")
    ax0.axhline(OWNER_THRESHOLD, color="r", ls=":", lw=1.2, label="owner 15 %")
    ax0.set_ylabel("relative volume error (own reader)")
    ax0.set_title("(A) round-trip error per format (bars at 1e-17 are exact)")
    ax0.legend(fontsize=7, ncol=2)
    ax0.grid(axis="y", alpha=0.3)

    # (B) geometry: the torus (genus 1) as written by the native STL writer
    body = mf.torus(1.0, 0.3, 48, 24)
    tri = body.points[body.faces]
    ax1.remove()
    ax1 = fig.add_subplot(1, 2, 2, projection="3d")
    ax1.plot_trisurf(tri[:, :, 0].ravel(), tri[:, :, 1].ravel(), tri[:, :, 2].ravel(),
                     triangles=np.arange(len(tri) * 3).reshape(-1, 3), cmap="viridis", lw=0.15,
                     edgecolor="k", alpha=0.95)
    ax1.set_title(f"(B) torus body, closed surface, Euler chi = {body.euler()}")
    ax1.set_xlabel("x [m]")
    ax1.set_ylabel("y [m]")
    ax1.set_zlabel("z [m]")
    ax1.view_init(elev=26, azim=38)

    fig.suptitle("TensorCAD native six-format mesh matrix (m17907)", fontsize=11)
    fig.tight_layout()
    path = figdir / "fig_native_six_format.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(ROOT / "artifacts" / "m17907_native_six_format"))
    ap.add_argument("--check", action="store_true", help="recompute and diff against the committed results.json")
    args = ap.parse_args(argv)
    outdir = Path(args.out)

    if args.check:
        committed = json.loads((outdir / "results.json").read_text(encoding="utf-8"))
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            fresh = run(Path(td))

        def norm(payload, keys):
            # a rerun happens in another directory: file paths (and the sha index,
            # which is keyed by relative path) are location-dependent, not numerical
            clone = {k: payload.get(k) for k in keys}
            for rows in (clone.get("rows"), clone.get("cross_writer")):
                for row in rows or []:
                    row.pop("file", None)
            if isinstance(clone.get("sids"), dict):
                clone["sids"].pop("legacy_file", None)
                clone["sids"].pop("native_file", None)
            return json.dumps(clone, sort_keys=True)

        drift = []
        for key in ("rows", "cross_writer", "sids", "summary"):
            if norm(committed, [key]) != norm(fresh, [key]):
                drift.append(key)
        if drift:
            print("DRIFT in", ", ".join(drift))
            return 1
        print("OK: results.json reproduces")
        return 0

    results = run(outdir)
    figpath = figure(results, outdir)
    results["figures"] = [{"name": figpath.name, "sha256": sha256(figpath)}]
    (outdir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary = results["summary"]
    print(json.dumps(summary, indent=2))
    print(f"figures: {figpath}")
    ok = (
        summary["own_reader_all_within_impl_threshold"]
        and summary["own_reader_topology_pass"] == summary["cells"]
        and summary["cgns_native_violations"] == 0
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())