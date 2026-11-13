#!/usr/bin/env python3
"""Cross-check of mission m17250 (anchor GAP-CAD-T1-IGES-PIPE) — audit mission m17306.

Non-destructive: reads the m17250 work tree (default /root/tp-agent/work/m17250) and
re-runs its published commands in a copy (/tmp/audit_m17250).  Every number in the
rendered report comes from ``checks.json`` produced here — nothing is hand-typed.

Usage: python tools/audit_m17306_crosscheck.py [src_work_tree]
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                      # /root/tp-agent/work/m17306
OUT = ROOT / "artifacts" / "audit_m17306"
SRC = Path(sys.argv[1] if len(sys.argv) > 1 else "/root/tp-agent/work/m17250").resolve()
ECO = Path("/root/tensor-ecosystem/TensorCAD")
DSN = ("host=127.0.0.1 port=5432 dbname=tensorplatform user=foundation "
       "password=uGTN-KwuQv_7l8Tw0MMWUx_kf-nF2Rzz")

DELIVERABLES = [
    "tensorcad/geomquality/iges_trim.py",
    "tensorcad/geomquality/trim_mesh.py",
    "tensorcad/geomquality/import_cli.py",
    "tests/test_iges_trim_import.py",
    "docs/gap_cad_t1_iges_pipe.md",
    "runs/gap_t1_iges/run_evidence.py",
    "runs/gap_t1_iges/manifest.json",
    "runs/gap_t1_iges/kvlcc2_class.report.json",
    "runs/gap_t1_iges/series60.report.json",
    "runs/gap_t1_iges/KP505_modelScale_250mm.report.json",
    "runs/gap_t1_iges/kvlcc2_class.stl",
    "runs/gap_t1_iges/series60.stl",
    "runs/gap_t1_iges/KP505_modelScale_250mm.stl",
    "artifacts/gap_t1_iges/provenance.json",
    "artifacts/gap_t1_iges/figure.png",
    "artifacts/gap_t1_iges/inputs/KP505_modelScale_250mm.igs",
]
REPORTS = {
    "KP505_modelScale_250mm.igs": "runs/gap_t1_iges/KP505_modelScale_250mm.report.json",
    "kvlcc2_class.iges": "runs/gap_t1_iges/kvlcc2_class.report.json",
    "series60.iges": "runs/gap_t1_iges/series60.report.json",
}
BASELINES = {"series60.iges": "series60", "kvlcc2_class.iges": "kvlcc2_class"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True).stdout


def inventory() -> list[dict]:
    head = set(git("ls-tree", "-r", "--name-only", "HEAD").splitlines())
    tracked_any = {}
    for rel in DELIVERABLES:
        log = git("log", "--all", "--oneline", "--", rel).strip()
        tracked_any[rel] = log.splitlines()[0] if log else None
    rows = []
    for rel in DELIVERABLES:
        p = SRC / rel
        rows.append({
            "path": rel,
            "exists": p.is_file(),
            "bytes": p.stat().st_size if p.is_file() else 0,
            "sha256": sha256(p) if p.is_file() else None,
            "in_main_HEAD": rel in head,
            "committed_anywhere": tracked_any[rel],
        })
    return rows


def benchmark_rows() -> dict:
    try:
        import psycopg

        with psycopg.connect(DSN) as conn, conn.cursor() as cur:
            cur.execute("SELECT count(*) FROM tensor_benchmark_cases")
            total = cur.fetchone()[0]
            cur.execute("SELECT case_name, metric, tolerance, repo, status FROM tensor_benchmark_cases "
                        "WHERE anchor_slug = 'GAP-CAD-T1-IGES-PIPE'")
            rows = cur.fetchall()
            cur.execute("SELECT anchor_slug, case_name, repo, status FROM tensor_benchmark_cases "
                        "WHERE anchor_slug ILIKE '%CAD%' OR case_name ILIKE '%IGES%' "
                        "OR case_name ILIKE '%STEP%' OR repo = 'TensorCAD'")
            others = cur.fetchall()
        return {"ok": True, "table_rows": total, "anchor_rows": rows,
                "related_rows": others, "anchor": "GAP-CAD-T1-IGES-PIPE"}
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "error": f"{type(exc).__name__}: {exc}"}


def reports() -> dict:
    out = {}
    for name, rel in REPORTS.items():
        d = json.loads((SRC / rel).read_text())
        out[name] = {
            "report": rel,
            "route": d.get("route", "trimmed"),
            "faces_in_file": d["import"].get("faces_in_file"),
            "faces_tessellated": d["import"].get("faces_tessellated"),
            "faces_refused": d["import"].get("faces_refused"),
            "refused_detail": d["import"].get("refused"),
            "max_inversion_residual": d["import"].get("max_inversion_residual"),
            "residual_tolerance": d["import"].get("residual_tolerance"),
            "entity_census": d["import"].get("entity_census"),
            "triangles": d["after_repair"]["triangles"],
            "free_edges_after": d["after_repair"]["free_edges"],
            "signed_volume": d["after_repair"]["signed_volume"],
            "watertight": d["after_repair"]["watertight"],
            "gates": {k: v["passed"] for k, v in d["gates"].items()},
            "gates_passed": all(v["passed"] for v in d["gates"].values()),
            "elapsed_seconds": d.get("elapsed_seconds"),
            "route_notes": d.get("route_notes", []),
        }
    manifest = json.loads((SRC / "runs/gap_t1_iges/manifest.json").read_text())
    out["_manifest_all_gates_passed"] = manifest.get("all_gates_passed")
    out["_manifest_files"] = [f["path"] for f in manifest["files"]]
    return out


def fidelity(rep: dict) -> dict:
    sys.path.insert(0, str(SRC))
    from tensorcad.geomquality import baselines  # noqa: PLC0415

    out = {}
    for name, fn in BASELINES.items():
        specimen = getattr(baselines, fn)()
        analytic = float(specimen.reference["volume"])
        mesh = float(rep[name]["signed_volume"])
        cb = float(specimen.particulars["block_coefficient"])
        out[name] = {
            "analytic_volume": analytic,
            "mesh_volume": mesh,
            "rel_deviation_pct": 100.0 * (mesh - analytic) / analytic,
            "analytic_cb": cb,
            "mesh_cb": cb * mesh / analytic,
            "length": specimen.length, "beam": specimen.beam, "draft": specimen.draft,
        }
    return out


def repro() -> dict:
    out = {"copy": "/tmp/audit_m17250"}
    march = "/tmp/audit_m17250/runs/repro_m17306"
    for name in ("series60.iges", "kvlcc2_class.iges"):
        stem = name[:-5]
        p = Path(f"{march}/{stem}.report.json")
        if p.is_file():
            d = json.loads(p.read_text())
            out[name] = {"volume": d["after_repair"]["signed_volume"],
                         "triangles": d["after_repair"]["triangles"]}
    kp = Path("/tmp/audit_m17250/runs/repro_kp505/KP505_modelScale_250mm.report.json")
    if kp.is_file():
        d = json.loads(kp.read_text())
        out["KP505_modelScale_250mm.igs"] = {
            "route": d.get("route", "trimmed"),
            "faces_tessellated": d["import"].get("faces_tessellated"),
            "faces_refused": d["import"].get("faces_refused"),
            "triangles": d["after_repair"]["triangles"],
            "refused_detail": d["import"].get("refused"),
        }
    return out


def tests() -> dict:
    log = OUT / "logs" / "pytest_iges.log"
    text = log.read_text() if log.is_file() else ""
    line = next((l for l in reversed(text.splitlines()) if "passed" in l), "")
    skip = next((l for l in text.splitlines() if l.startswith("SKIPPED")), "")
    return {"log": str(log.relative_to(ROOT)), "summary": line,
            "skipped": skip, "skipif_source": "tests/test_iges_trim_import.py:203"}


def gold() -> dict:
    p = OUT / "gold_audit.json"
    if not p.is_file():
        return {"status": "not_run"}
    d = json.loads(p.read_text())
    out = {"status": "run"}
    for key in ("gold_exact", "plain"):
        e = d.get(key, {})
        if "error" in e:
            out[key] = e
            continue
        out[key] = {"raw_faces": e["raw_faces"], "raw_vertices": e["raw_vertices"],
                    "vertices_merged": e["vertices_merged"],
                    "free_edges_before": e["before"]["free_edges"],
                    "watertight_before": e["before"]["watertight"],
                    "signed_volume_before": e["before"]["signed_volume"],
                    "free_edges_after_repair": e["after_repair"]["free_edges"],
                    "watertight_after_repair": e["after_repair"]["watertight"],
                    "signed_volume_after_repair": e["after_repair"]["signed_volume"],
                    "trimesh_is_watertight": e["trimesh_is_watertight"]}
    return out


def integration() -> dict:
    """Does anything in the shipped package import the new modules?"""
    hits = subprocess.run(
        ["grep", "-rn", "--include=*.py", "-E",
         r"(geomquality\.)?(trim_mesh|iges_trim|import_cli)", str(ECO / "tensorcad")],
        capture_output=True, text=True).stdout.splitlines()
    hits = [h for h in hits if "/geomquality/trim_mesh.py" not in h
            and "/geomquality/iges_trim.py" not in h and "/geomquality/import_cli.py" not in h]
    return {"eco_repo": str(ECO), "importers_outside_new_modules": hits[:10],
            "importer_count": len(hits)}


def main() -> dict:
    rep = reports()
    checks = {
        "anchor": "GAP-CAD-T1-IGES-PIPE",
        "mission_under_audit": 17250,
        "audit_mission": 17306,
        "src_work_tree": str(SRC),
        "inventory": inventory(),
        "benchmark": benchmark_rows(),
        "reports": rep,
        "fidelity": fidelity(rep),
        "repro": repro(),
        "tests": tests(),
        "gold": gold(),
        "integration": integration(),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "checks.json").write_text(json.dumps(checks, indent=2, default=str))
    print(json.dumps({k: checks[k] for k in ("benchmark", "fidelity", "integration")}, indent=1)[:2000])
    return checks


if __name__ == "__main__":
    main()