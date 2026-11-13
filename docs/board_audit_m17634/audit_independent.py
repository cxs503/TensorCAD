#!/usr/bin/env python3
"""Independent audit recompute for m17595 (XVAL-CAD-M17483-NATIVE-SIX-FORMAT-MATRIX).

Independence rule: this script MUST NOT import ``tensorcad.mesh_xval`` nor the
subject's tool. Geometry is re-read with third-party readers only (meshio 5.3.5,
trimesh, numpy-stl, h5py), plus a ~25-line SU2 reader written here from the SU2
v7 grammar, so every number below is a second opinion, not a re-run.

Usage: python3 docs/board_audit_m17634/audit_independent.py [subject_root]
       (default subject_root = /tmp/m17595_verify, a git-archive-free clone of
        branch agent/m17595; see docs/board_audit_m17634/audit_report.md)

Outputs: docs/board_audit_m17634/audit_evidence.json (+ prints a summary).
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_SUBJECT = Path("/tmp/m17595_verify")
ART_REL = "artifacts/m17595_xval_six_format"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def soup_volume(tris: np.ndarray) -> float:
    """Signed volume of a triangle soup via the divergence theorem."""
    v0, v1, v2 = tris[:, 0], tris[:, 1], tris[:, 2]
    return float(abs(np.sum(np.einsum("ij,ij->i", v0, np.cross(v1, v2))) / 6.0))


def read_su2(path: Path) -> np.ndarray:
    """Minimal SU2 mesh-file reader written for this audit.

    Collects triangles from the volume section (``NELEM``) and from every boundary
    marker block (``MARKER_ELEMS``), so the ``NELEM=0`` surface-mesh dialect that the
    subject ships as a variant is read the same way a solver's surface wrapper would.
    """
    toks = [ln.split("%")[0].strip() for ln in Path(path).read_text().splitlines()]
    toks = [ln for ln in toks if ln]
    pos = {}
    for k, ln in enumerate(toks):
        if "=" in ln:
            pos.setdefault(ln.split("=")[0].strip(), k)
    npoin = int(toks[pos["NPOIN"]].split("=")[1])
    pts = np.array([[float(x) for x in toks[pos["NPOIN"] + 1 + k].split()[:3]]
                    for k in range(npoin)])
    tris = []
    for key in ("NELEM", "MARKER_ELEMS"):
        k = pos.get(key)
        if k is None:
            continue
        n = int(toks[k].split("=")[1])
        for j in range(n):
            parts = toks[k + 1 + j].split()
            if int(parts[0]) != 5:  # 5 == triangle in the SU2 element type table
                raise ValueError(f"{path.name}: {key} entry {j} is type {parts[0]}, not a triangle")
            tris.append([int(parts[1]), int(parts[2]), int(parts[3])])
    return pts[np.array(tris, dtype=int)]


def read_cgns(label_based: bool):
    """Audit-side CGNS readers, written here from the SIDS-to-HDF5 mapping.

    ``label_based=True`` follows the SIDS-to-HDF5 convention (node identity in the
    ``name``/``label`` attributes, node names free); ``False`` follows the
    name-addressed dialect that the upstream m17483 files (and meshio's writer)
    use.  Neither path touches ``tensorcad.mesh_xval``.
    """
    import h5py

    def _attr(g, key):
        for k in (key, " " + key):
            if k in g.attrs:
                return g.attrs[k]
        return b""

    def _label(g):
        v = _attr(g, "label")
        return v.decode() if isinstance(v, bytes) else str(v)

    def _name(g, key):
        v = _attr(g, "name")
        v = v.decode() if isinstance(v, bytes) else str(v)
        return v or key

    def read(path: Path) -> np.ndarray:
        with h5py.File(path, "r") as f:
            base = next((f[k] for k in f.keys() if isinstance(f[k], h5py.Group)
                         and (not label_based or _label(f[k]) == "CGNSBase_t")), None)
            if base is None:
                raise ValueError("no CGNSBase")
            if label_based:
                zone = next(g for k in base.keys() if isinstance((g := base[k]), h5py.Group)
                            and _label(g) == "Zone_t")
                coords = next(zone[k] for k in zone.keys() if isinstance(zone[k], h5py.Group)
                              and _label(zone[k]) == "GridCoordinates_t")
                elems = [zone[k] for k in zone.keys() if isinstance(zone[k], h5py.Group)
                         and _label(zone[k]) == "Elements_t"]
            else:
                zone = base["Zone1"]
                coords = zone["GridCoordinates"]
                elems = [zone["GridElements"]]
            def arr(node):
                d = node[" data"]
                return np.asarray(d, dtype=float)
            xyz = [arr(coords[c]) for c in ("CoordinateX", "CoordinateY", "CoordinateZ")]
            pts = np.column_stack(xyz)
            tris = []
            for e in elems:
                etype = int(np.asarray(e["data" if "data" in e else " data"]).ravel()[0])
                conn = np.asarray(e["ElementConnectivity"][" data"]).ravel()
                lo, hi = np.asarray(e["ElementRange"][" data"]).ravel()
                n = int(hi - lo + 1)
                if etype == 5:      # triangle
                    tris.append(conn[: 3 * n].reshape(n, 3))
                elif etype == 7:    # quad -> split into two triangles
                    q = conn[: 4 * n].reshape(n, 4)
                    tris.append(q[:, [0, 1, 2]])
                    tris.append(q[:, [0, 2, 3]])
                else:
                    raise ValueError(f"element type {etype} unsupported")
            return pts[(np.vstack(tris) - 1)]

    return read

def read_with(lib: str, path: Path) -> np.ndarray:
    if lib == "meshio":
        import meshio

        m = meshio.read(path)
        blocks = [c for c in m.cells if c.type == "triangle"]
        if not blocks:
            raise ValueError("no triangle block")
        idx = np.vstack([b.data for b in blocks])
        pts = np.asarray(m.points, dtype=float)
        return pts[idx]
    if lib == "trimesh":
        import trimesh

        m = trimesh.load(str(path), force="mesh", process=False)
        return np.asarray(m.vertices)[np.asarray(m.faces)]
    if lib == "numpy-stl":
        from stl import mesh as stlmesh

        m = stlmesh.Mesh.from_file(str(path))
        return np.asarray(m.vectors, dtype=float)
    if lib == "su2-audit":
        return read_su2(path)
    if lib == "cgns-sids-audit":
        return read_cgns(True)(path)
    if lib == "cgns-namebased-audit":
        return read_cgns(False)(path)
    raise ValueError(lib)


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SUBJECT)
    art = root / ART_REL
    res = json.loads((art / "results.json").read_text())
    out: dict = {"subject_root": str(root), "checks": {}}

    # ── C1: artefact manifest integrity (sha256 of every shipped file) ────────
    man = json.loads((art / "manifest.json").read_text())
    missing, mismatch, checked = [], [], 0
    for rel, meta in man["files"].items():
        p = art / rel
        if not p.is_file():
            missing.append(rel)
            continue
        checked += 1
        if sha256(p) != meta["sha256"]:
            mismatch.append(rel)
    doc = root / man["doc"]["path"]
    doc_ok = doc.is_file() and sha256(doc) == man["doc"]["sha256"]
    out["checks"]["manifest"] = {
        "files_declared": len(man["files"]), "verified": checked,
        "missing": missing[:5], "sha256_mismatch": mismatch[:5],
        "doc_path": man["doc"]["path"], "doc_sha256_matches": doc_ok,
    }

    # ── C2: independent re-read of the six TensorCAD-written formats ──────────
    # reader set per format is chosen by the third-party availability only.
    readers = {
        "stl": ["meshio", "trimesh", "numpy-stl"],
        "obj": ["meshio", "trimesh"],
        "vtk": ["meshio"],
        "msh": ["meshio"],
        "cgns": ["meshio", "cgns-sids-audit"],
        "su2": ["su2-audit"],
    }
    cells = {}
    for gdir in sorted((art / "files").iterdir()):
        if not gdir.is_dir():
            continue
        for fmt, libs in readers.items():
            p = gdir / f"{gdir.name}.{fmt}"
            if not p.is_file():
                continue
            for lib in libs:
                key = f"{gdir.name}/{fmt}/{lib}"
                try:
                    tris = read_with(lib, p)
                    cells[key] = {"ok": True, "volume": soup_volume(tris),
                                  "triangles": int(len(tris)), "bytes": p.stat().st_size}
                except BaseException as e:  # unreadable by that library = honest finding
                    cells[key] = {"ok": False, "error": f"{type(e).__name__}: {str(e)[:160]}",
                                  "bytes": p.stat().st_size}
    out["checks"]["independent_reread"] = cells

    # ── C3: absolute + cross-format agreement per geometry ───────────────────
    analytic = {"cube": 8.0}
    agree = {}
    for gdir in sorted((art / "files").iterdir()):
        if not gdir.is_dir():
            continue
        vols = {k: v["volume"] for k, v in cells.items() if k.startswith(gdir.name + "/") and v["ok"]}
        if not vols:
            continue
        vals = np.array(list(vols.values()))
        ref_v = analytic.get(gdir.name, float(np.median(vals)))
        agree[gdir.name] = {
            "formats_read_independently": sorted({k.split("/")[1] for k in vols}),
            "readers": len(vols),
            "volumes": {k: round(v, 12) for k, v in vols.items()},
            "reference": ref_v,
            "max_rel_err_vs_reference": float(np.max(np.abs(vals - ref_v)) / ref_v),
            "cross_format_spread_rel": float((vals.max() - vals.min()) / max(vals.min(), 1e-30)),
            "analytic_exact": gdir.name in analytic,
        }
    out["checks"]["volume_agreement"] = agree

    # ── C3b: independent re-read of the upstream m17483 files (matrix B) ─────
    up = art / "upstream_m17483" / "files"
    upstream = {}
    for p in sorted(up.iterdir()):
        fmt = p.suffix.lstrip(".")
        if fmt in ("stl", "obj", "vtk", "msh"):
            try:
                upstream[p.name] = {"ok": True, "volume": soup_volume(read_with("meshio", p))}
            except BaseException as e:
                upstream[p.name] = {"ok": False, "error": f"{type(e).__name__}: {str(e)[:120]}"}
        elif fmt == "cgns":
            for lib in ("cgns-namebased-audit", "cgns-sids-audit", "meshio"):
                try:
                    upstream[p.name + "#" + lib] = {"ok": True,
                                                   "volume": soup_volume(read_with(lib, p))}
                except BaseException as e:
                    upstream[p.name + "#" + lib] = {"ok": False,
                                                    "error": f"{type(e).__name__}: {str(e)[:120]}"}
        elif fmt == "su2":
            try:
                upstream[p.name] = {"ok": True, "volume": soup_volume(read_with("su2-audit", p))}
            except BaseException as e:
                upstream[p.name] = {"ok": False, "error": f"{type(e).__name__}: {str(e)[:120]}"}
    out["checks"]["upstream_reread"] = upstream

    # ── C4: does the subject's headline number survive recount? ──────────────
    claimed, status_hist = [], {}
    for mkey in ("matrix_a", "matrix_b"):
        m = res.get(mkey, {})
        for cell in m.get("cells", []):
            for rd in cell.get("readers", []):
                st = rd.get("status", "?")
                status_hist[st] = status_hist.get(st, 0) + 1
                claimed.append({"matrix": mkey, "geometry": cell.get("geometry"),
                                "format": cell.get("format"), "reader": rd.get("reader"),
                                "variant": rd.get("variant"), "status": st,
                                "volume_rel_err": rd.get("volume_rel_err")})
    errs = [c["volume_rel_err"] for c in claimed if c["volume_rel_err"] is not None]
    out["checks"]["claimed_errors"] = {
        "reader_checks_total": len(claimed),
        "status_histogram": status_hist,
        "checks_with_a_number": len(errs),
        "max_volume_rel_err": max(errs) if errs else None,
        "worst_cell": max((c for c in claimed if c["volume_rel_err"] is not None),
                          key=lambda c: c["volume_rel_err"]) if errs else None,
        "owner_threshold": res["acceptance"]["owner_threshold"],
        "implementation_threshold": res["acceptance"]["implementation_threshold"],
        "all_numbered_checks_within_impl_threshold":
            all(e <= res["acceptance"]["implementation_threshold"] for e in errs),
        "headline_claim": "worst 1.311e-8 (18 cells, 33+30 reader checks)",
        "headline_reproduced": bool(errs) and abs(max(errs) - 1.311e-08) / 1.311e-08 < 1e-3,
        "summary_A": res.get("matrix_a", {}).get("summary", {}),
        "summary_B_declared_reader_checks": res.get("matrix_b", {}).get("summary", {}).get(
            "reader_checks"),
        "summary_B_declared_green": res.get("matrix_b", {}).get("summary", {}).get(
            "reader_checks_green"),
    }

    # ── C5: CGNS — verify the upstream-mapping defect claim with h5py ─────────
    import h5py

    def cgns_probe(p: Path) -> dict:
        """Count SIDS-to-HDF5 mapping violations, using the standard's own rules:
        every node carries name/label/order attributes and every data node stores
        its payload in a dataset called " data"; a CGNSLibraryVersion node with a
        " version" attribute must sit at the root."""
        info = {"file": str(p.relative_to(root)), "nodes": 0, "missing_name": 0,
                "missing_label": 0, "missing_order": 0, "payload_as_dataset": 0,
                "label_in_type_attr": 0, "nonstandard_order": 0,
                "root_has_cgnslibversion_node": False, "root_has_version_attr": False}

        def attrs(obj, key):
            for k in (key, " " + key):
                if k in obj.attrs:
                    v = obj.attrs[k]
                    return v.decode() if isinstance(v, bytes) else v
            return None

        with h5py.File(p, "r") as f:
            info["root_has_version_attr"] = attrs(f, "version") not in (None, "")
            for k in f.keys():
                info["root_has_cgnslibversion_node"] |= k.strip() == "CGNSLibraryVersion"

            def walk(g):
                for k in g.keys():
                    o = g[k]
                    if not isinstance(o, h5py.Group):
                        info.setdefault("payload_datasets", 0)
                        info["payload_datasets"] += 1
                        continue
                    info["nodes"] += 1
                    if attrs(o, "name") in (None, ""):
                        info["missing_name"] += 1
                    if attrs(o, "label") in (None, ""):
                        info["missing_label"] += 1
                        if attrs(o, "type") not in (None, ""):
                            info["label_in_type_attr"] += 1
                    if attrs(o, "order") in (None, ""):
                        info["missing_order"] += 1
                    if "data" in o.keys() or " data" in o.keys():
                        info["payload_as_dataset"] += 1
                    walk(o)

            walk(f)
        return info

    probe = {}
    for tag, p in (("xval_written_cube", art / "files" / "cube" / "cube.cgns"),
                   ("upstream_m17483_cube", art / "upstream_m17483" / "files" / "cube.cgns")):
        if p.is_file():
            try:
                probe[tag] = cgns_probe(p)
            except Exception as e:
                probe[tag] = {"error": f"{type(e).__name__}: {str(e)[:120]}"}
    out["checks"]["cgns_probe"] = probe

    # ── C6: report six-section completeness + figure inventory ───────────────
    doc_txt = doc.read_text() if doc.is_file() else ""
    out["checks"]["report"] = {
        "doc": man["doc"]["path"],
        "sections": [ln.strip() for ln in doc_txt.splitlines()
                     if ln.startswith("## ") and not ln.startswith("### ")],
        "has_error_table": "15 %" in doc_txt or "15%" in doc_txt,
        "declares_nonzero_mismatch": "NO — see" in doc_txt,
        "figures": {p.name: p.stat().st_size for p in sorted((art / "figures").glob("*.png"))},
    }

    (HERE / "audit_evidence.json").write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: (v if k != "independent_reread" else f"{len(v)} reader entries")
                      for k, v in out["checks"].items()}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())