"""GAP-CAD-T1-IGES-PIPE benchmark: IGES/STEP -> watertight mesh import pipeline.

Reproduce (from the repository root)::

    python benchmarks/verified/iges_import/run_benchmark.py --out-dir runs/gap_t1_iges

What it does
------------
1. writes the input CAD files (``artifacts/gap_t1_iges/inputs``) from the repository's
   analytic hull baselines and from the platform's ``DTMB5415_Hull_exact.stl`` gold mesh
   (``gold_loft.py``), recording provenance (path, sha256, size) for each;
2. imports every file through ``tensorcad.geometry.load_cad_geometry`` -- the public
   integration entry point, which calls ``geomquality.trim_mesh.import_cad`` -- and runs
   the audit/repair/gate chain;
3. benchmarks *ref* against *got*: the reference volume is a 5x5 Gauss-Legendre
   quadrature over the very NURBS patches that were written (independent of
   triangulation), the measured volume is the repaired triangle mesh that the pipeline
   returns.  Deviation of the DTMB5415 round trip against the gold STL is measured on
   transverse sections (naval-architecture convention);
4. sweeps the patch-set tessellation density to show the discretisation error is a
   closable ``O(h^2)`` term rather than a fixed bias.

Nothing in this file fabricates a number: every entry of ``evidence.json`` is read back
from the run it just performed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from gold_loft import (  # noqa: E402
    brep_quadrature, brep_quadrature as _quad, loft_gold_hull, read_stl, section_polyline,
)

from tensorcad.geomquality.baselines import kvlcc2_class, series60  # noqa: E402
from tensorcad.geomquality.io import write_iges  # noqa: E402
from tensorcad.geomquality.trim_mesh import import_cad  # noqa: E402
from tensorcad.geometry import load_cad_geometry  # noqa: E402

GOLD = Path("/root/ship-performance-platform/backend/data/geometry/ships/DTMB5415_Hull_exact.stl")
KP505 = Path("/root/tensor-ecosystem/TensorBEM/benchmarks/anchors/sources/KP505_modelScale_250mm.igs")
TOLERANCE = 3.0  # per-cent, the acceptance threshold of the anchor


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _mesh_volume(vertices, faces) -> float:
    points = np.asarray(vertices, dtype=float)
    triangles = np.asarray(faces, dtype=int)
    a, b, c = points[triangles[:, 0]], points[triangles[:, 1]], points[triangles[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def _section_deviation(gold_triangles: np.ndarray, vertices, faces, *, stations: int = 15):
    """Mean/max distance between the imported and gold transverse sections, in metres.

    Both meshes are cut at the same stations; each section is resampled on arclength and
    cropped to ``z <= min(z_max_import, z_max_gold) - 0.02 m`` so the comparison covers
    the wetted girth and not the deck/cap closure faces.
    """
    points = np.asarray(vertices, dtype=float)
    triangles = np.asarray(faces, dtype=int)
    mesh_triangles = points[triangles]
    x_min = max(gold_triangles[:, :, 0].min(), mesh_triangles[:, :, 0].min())
    x_max = min(gold_triangles[:, :, 0].max(), mesh_triangles[:, :, 0].max())
    length = x_max - x_min
    distances: list[float] = []
    for fraction in np.linspace(0.05, 0.95, stations):
        x0 = x_min + fraction * length
        gold_chains = section_polyline(gold_triangles, float(x0))
        mesh_chains = section_polyline(mesh_triangles, float(x0))
        if not gold_chains or not mesh_chains:
            continue
        gold = max(gold_chains, key=lambda c: np.ptp(c[:, 1]))
        mesh = max(mesh_chains, key=lambda c: np.ptp(c[:, 1]))
        z_limit = min(gold[:, 2].max(), mesh[:, 2].max()) - 0.02
        gold = gold[gold[:, 2] <= z_limit]
        mesh = mesh[mesh[:, 2] <= z_limit]
        if len(gold) < 4 or len(mesh) < 4:
            continue
        # resample both on normalised arclength, then nearest-distance (y, z)
        gold_s = _resample(gold, 240)
        mesh_s = _resample(mesh, 240)
        delta = gold_s[:, 1:, None] - mesh_s[None, :, 1:]
        distances.extend(np.sqrt((delta ** 2).sum(axis=2)).min(axis=1).tolist())
    if not distances:
        return {"mean": None, "max": None, "pct_L": None, "stations": 0}
    values = np.asarray(distances)
    return {
        "mean": float(values.mean()),
        "max": float(values.max()),
        "pct_L": float(100.0 * values.mean() / (x_max - x_min)),
        "stations": stations,
        "note": "distance from gold sections to imported sections, wetted girth only",
    }


def _resample(points: np.ndarray, count: int) -> np.ndarray:
    steps = np.linalg.norm(np.diff(points, axis=0), axis=1)
    total = steps.sum()
    if total <= 0.0:
        return points
    cumulative = np.concatenate([[0.0], np.cumsum(steps)]) / total
    targets = np.linspace(0.0, 1.0, count)
    return np.stack([np.interp(targets, cumulative, points[:, axis]) for axis in range(3)], axis=1)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", default="runs/gap_t1_iges")
    parser.add_argument("--samples", default="12x8,24x12,48x24,96x48",
                        help="convergence sweep, <u>x<v> pairs")
    args = parser.parse_args(argv)

    out = ROOT / args.out_dir
    inputs = ROOT / "artifacts" / "gap_t1_iges" / "inputs"
    out.mkdir(parents=True, exist_ok=True)
    inputs.mkdir(parents=True, exist_ok=True)

    started = time.time()
    provenance: dict[str, dict] = {}
    brep_refs: dict[str, dict] = {}
    gold_triangles = None

    # ---- input files -------------------------------------------------------
    builders = {
        "series60_cb060": lambda: series60(0.60),
        "series60_cb065": lambda: series60(0.65),
        "kvlcc2_class": kvlcc2_class,
    }
    for name, builder in builders.items():
        specimen = builder()
        target = inputs / f"{name}.iges"
        write_iges(target, specimen.brep)
        provenance[name] = {
            "path": str(target.relative_to(ROOT)), "sha256": sha256(target),
            "bytes": target.stat().st_size, "kind": "analytic hull baseline",
            "source": "tensorcad.geomquality.baselines", "license": "in-repository",
        }
        brep_refs[name] = brep_quadrature(specimen.brep)

    if GOLD.exists():
        gold_triangles = read_stl(GOLD)
        brep, info = loft_gold_hull(gold_triangles)
        target = inputs / "DTMB5415_gold_loft.iges"
        write_iges(target, brep)
        provenance["DTMB5415_gold_loft"] = {
            "path": str(target.relative_to(ROOT)), "sha256": sha256(target),
            "bytes": target.stat().st_size, "kind": "NURBS loft of the platform gold mesh",
            "source": str(GOLD), "source_sha256": sha256(GOLD),
            "source_faces": int(len(gold_triangles)),
            "license": "platform internal reference geometry", "loft": info,
        }
        brep_refs["DTMB5415_gold_loft"] = brep_quadrature(brep)
    else:
        print(f"WARNING gold mesh missing: {GOLD}")

    # ---- import + measure --------------------------------------------------
    rows: dict[str, dict] = {}
    for name, reference in brep_refs.items():
        path = inputs / f"{name}.iges"
        record: dict = {"input": str(path.relative_to(ROOT)), "sha256": provenance[name]["sha256"]}
        try:
            clock = time.time()
            mesh, report = load_cad_geometry(path)
            record["seconds"] = round(time.time() - clock, 2)
            volume = _mesh_volume(mesh.vertices, mesh.faces)
            after = report["after_repair"]
            record.update({
                "route": report.get("route", "trimmed"),
                "triangles": int(after["triangles"]),
                "watertight": bool(after["watertight"]),
                "free_edges": int(after["free_edges"]),
                "components": int(after["component_count"]),
                "gates": {key: bool(value["passed"]) for key, value in report["gates"].items()},
                "gates_passed": all(value["passed"] for value in report["gates"].values()),
                "ref_volume_quadrature": reference["volume"],
                "got_volume_mesh": volume,
                "volume_error_pct": 100.0 * (volume - reference["volume"]) / reference["volume"],
                "ref_area_quadrature": reference["area"],
                "got_area_mesh": float(after.get("area", float("nan"))),
                "area_error_pct": (100.0 * (float(after.get("area", 0.0)) - reference["area"])
                                   / reference["area"]),
                "within_3pct": abs(100.0 * (volume - reference["volume"]) / reference["volume"]) < TOLERANCE,
            })
            if gold_triangles is not None and name == "DTMB5415_gold_loft":
                record["deviation_vs_gold"] = _section_deviation(gold_triangles, mesh.vertices, mesh.faces)
            (out / f"{name}.report.json").write_text(json.dumps(report, indent=2, default=str))
        except Exception as exc:  # honest failure, never a fabricated row
            record.update({"error": f"{type(exc).__name__}: {exc}", "gates_passed": False,
                           "within_3pct": False})
        rows[name] = record
        print(name, json.dumps({k: record.get(k) for k in
                                ("triangles", "watertight", "gates_passed", "volume_error_pct",
                                 "within_3pct")}))

    # ---- KP505: real third-party CATIA IGES (trimmed route) ----------------
    kp505: dict = {"path": str(KP505)}
    if KP505.exists():
        try:
            clock = time.time()
            result = import_cad(KP505, allow_patch_fallback=False)
            report = result.report
            kp505.update({
                "sha256": sha256(KP505), "seconds": round(time.time() - clock, 2),
                "faces_in_file": report["import"]["faces_in_file"],
                "faces_tessellated": report["import"]["faces_tessellated"],
                "faces_refused": report["import"]["faces_refused"],
                "max_inversion_residual": report["import"]["max_inversion_residual"],
                "refused": report["import"]["refused"],
                "triangles": int(report["after_repair"]["triangles"]),
                "watertight": bool(report["after_repair"]["watertight"]),
                "gates": {key: bool(value["passed"]) for key, value in report["gates"].items()},
            })
        except Exception as exc:
            kp505["error"] = f"{type(exc).__name__}: {exc}"
        print("KP505", json.dumps({k: kp505.get(k) for k in
                                  ("faces_in_file", "faces_tessellated", "faces_refused", "error")}))

    # ---- convergence sweep -------------------------------------------------
    sweep: dict[str, list] = {}
    pairs = []
    for token in args.samples.split(","):
        u, _, v = token.partition("x")
        pairs.append((int(u), int(v)))
    for name in ("series60_cb060", "DTMB5415_gold_loft"):
        if name not in brep_refs:
            continue
        path = inputs / f"{name}.iges"
        series = []
        for u, v in pairs:
            try:
                mesh, report = load_cad_geometry(path, u_samples=u, v_samples=v)
                volume = _mesh_volume(mesh.vertices, mesh.faces)
                ref = brep_refs[name]["volume"]
                series.append({
                    "u_samples": u, "v_samples": v,
                    "triangles": int(report["after_repair"]["triangles"]),
                    "volume": volume,
                    "error_pct": 100.0 * (volume - ref) / ref,
                    "watertight": bool(report["after_repair"]["watertight"]),
                })
            except Exception as exc:
                series.append({"u_samples": u, "v_samples": v, "error": str(exc)})
        sweep[name] = series
        print("convergence", name, [(s["u_samples"], s["v_samples"],
                                     round(s.get("error_pct", float("nan")), 4)) for s in series])

    evidence = {
        "anchor": "GAP-CAD-T1-IGES-PIPE",
        "tool": "tensorcad.geometry.load_cad_geometry -> geomquality.trim_mesh.import_cad",
        "tolerance_pct": TOLERANCE,
        "reference": "5x5 Gauss-Legendre quadrature over the written NURBS patches",
        "gold_mesh": ({"path": str(GOLD), "sha256": sha256(GOLD), "faces": int(len(gold_triangles))}
                      if gold_triangles is not None else None),
        "provenance": provenance,
        "files": rows,
        "kp505_real_catia": kp505,
        "convergence": sweep,
        "elapsed_seconds": round(time.time() - started, 1),
        "repro_cmd": "python benchmarks/verified/iges_import/run_benchmark.py --out-dir runs/gap_t1_iges",
    }
    (out / "evidence.json").write_text(json.dumps(evidence, indent=2, default=str))
    (ROOT / "artifacts" / "gap_t1_iges" / "provenance.json").write_text(
        json.dumps({"inputs": provenance, "gold_mesh": evidence["gold_mesh"]}, indent=2, default=str)
    )
    print(json.dumps({"evidence": str(out / "evidence.json"),
                      "all_gates_passed": all(r.get("gates_passed") for r in rows.values()),
                      "all_within_3pct": all(r.get("within_3pct") for r in rows.values())}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())