#!/usr/bin/env python3
"""Measure / publish stages for audit mission m17306 (m17250 cross-check).

stage measure : quick independent re-measurement of the two DTMB5415 STL files quoted in
                doc section 5 (fast: trimesh load + face counts + enclosed volume), and
                the audit figure.
stage publish : INSERT the maturity verdict into tensor_maturity_verdicts and POST the
                paper-grade report to the platform.
"""
from __future__ import annotations

import base64
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path("/root/tp-agent/work/m17306")
OUT = ROOT / "artifacts" / "audit_m17306"
DOCS = ROOT / "docs" / "audit_m17306"
DSN = ("host=127.0.0.1 port=5432 dbname=tensorplatform user=foundation "
       "password=uGTN-KwuQv_7l8Tw0MMWUx_kf-nF2Rzz")
API = "http://127.0.0.1:8010/api/v1"
ENV = dict(l.strip().split("=", 1) for l in open("/root/tp-agent/agent.env")
           if l.strip() and not l.startswith("#") and "=" in l)
STLS = {
    "gold_exact": "/root/ship-performance-platform/backend/data/geometry/ships/DTMB5415_Hull_exact.stl",
    "plain": "/root/ship-performance-platform/backend/data/geometry/ships/DTMB5415_Hull.stl",
}


def measure() -> None:
    import trimesh

    quick = {}
    for key, path in STLS.items():
        mesh = trimesh.load(path, process=False)
        quick[key] = {
            "path": path,
            "faces": int(len(mesh.faces)),
            "vertices": int(len(mesh.vertices)),
            "trimesh_volume": float(mesh.volume),
            "trimesh_is_watertight": bool(mesh.is_watertight),
            "trimesh_is_winding_consistent": bool(mesh.is_winding_consistent),
            "extent": [float(x) for x in mesh.extents],
        }
    (OUT / "gold_quick.json").write_text(json.dumps(quick, indent=2))

    checks = json.loads((OUT / "checks.json").read_text())
    quick["_gold_audit_slow"] = checks.get("gold", {}).get("status")
    out = quick
    (OUT / "gold_quick.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=1))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.0))
    names = list(checks["fidelity"])
    dev = [checks["fidelity"][n]["rel_deviation_pct"] for n in names]
    cb_a = [checks["fidelity"][n]["analytic_cb"] for n in names]
    cb_m = [checks["fidelity"][n]["mesh_cb"] for n in names]
    axes[0].bar(names, dev, color="#c44e52")
    axes[0].axhline(-3.0, ls="--", lw=1, color="k")
    axes[0].axhline(0.0, lw=0.8, color="k")
    axes[0].set_ylabel("imported mesh volume vs analytic hull [%]")
    axes[0].set_title("volume fidelity of the IGES round trip")
    for i, v in enumerate(dev):
        axes[0].text(i, v - 0.15, f"{v:.2f}%", ha="center", va="top", fontsize=9)
    width = 0.35
    idx = range(len(names))
    axes[1].bar([i - width / 2 for i in idx], cb_a, width, label="analytic Cb (baseline)")
    axes[1].bar([i + width / 2 for i in idx], cb_m, width, label="Cb of imported mesh")
    axes[1].set_xticks(list(idx))
    axes[1].set_xticklabels(names, fontsize=9)
    axes[1].set_ylabel("Cb")
    axes[1].set_title("Cb from the watertight gate output")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    png = OUT / "figure.png"
    fig.savefig(png, dpi=140)
    print("figure", png)


def verdict_evidence() -> dict:
    checks = json.loads((OUT / "checks.json").read_text())
    rep = checks["reports"]
    fid = checks["fidelity"]
    kp = rep["KP505_modelScale_250mm.igs"]
    return {
        "impl_files": [
            "/root/tp-agent/work/m17250/tensorcad/geomquality/iges_trim.py:1 (46871 B)",
            "/root/tp-agent/work/m17250/tensorcad/geomquality/trim_mesh.py:1 (43814 B)",
            "/root/tp-agent/work/m17250/tensorcad/geomquality/import_cli.py:1 (6524 B)",
        ],
        "entry": "/root/tp-agent/work/m17250/tensorcad/geomquality/import_cli.py (python -m "
                 "tensorcad.geomquality.import_cli <file.igs> --out-dir ... --allow-open)",
        "run_artifacts": [
            "/root/tp-agent/work/m17250/runs/gap_t1_iges/kvlcc2_class.report.json",
            "/root/tp-agent/work/m17250/runs/gap_t1_iges/series60.report.json",
            "/root/tp-agent/work/m17250/runs/gap_t1_iges/KP505_modelScale_250mm.report.json",
            "/root/tp-agent/work/m17250/runs/gap_t1_iges/kvlcc2_class.stl (451911 B)",
            "/root/tp-agent/work/m17250/runs/gap_t1_iges/KP505_modelScale_250mm.stl (49 B = empty, 0 triangles)",
        ],
        "tests": [
            "tests/test_iges_trim_import.py::test_plane_inversion_is_exact",
            "tests/test_iges_trim_import.py::test_revolved_surface_inversion_recovers_the_generating_parameters",
            "tests/test_iges_trim_import.py::test_triangulation_of_a_square_has_exact_area",
            "tests/test_iges_trim_import.py::test_triangulation_with_a_hole_removes_the_hole_area",
            "tests/test_iges_trim_import.py::test_ear_clip_keeps_the_area_of_a_non_convex_loop",
            "tests/test_iges_trim_import.py::test_closed_box_audits_exactly",
            "tests/test_iges_trim_import.py::test_repair_chain_closes_a_box_with_one_face_missing",
            "tests/test_iges_trim_import.py::test_weld_merges_duplicate_corner_vertices_only_within_tolerance",
            "run: 10 passed, 1 skipped (tests/test_iges_trim_import.py:203 needs runs/gap_t1_iges/evidence.json, which is absent), 1 xfailed (hole bridging)",
        ],
        "integration_callers": None,
        "integration_note": "grep -rn 'trim_mesh|iges_trim|import_cli' over /root/tensor-ecosystem/TensorCAD/tensorcad -> 0 importers outside the new modules'",
        "benchmark": None,
        "benchmark_note": "SELECT count(*) FROM tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-T1-IGES-PIPE' -> 0 rows (table has 46 rows); no ref/got/tolerance exists, so no error statement is possible",
        "repro_cmd": "cp -a /root/tp-agent/work/m17250 /tmp/audit_m17250 && cd /tmp/audit_m17250 && "
                     "python -m tensorcad.geomquality.import_cli artifacts/gap_t1_iges/inputs/"
                     "series60.iges artifacts/gap_t1_iges/inputs/kvlcc2_class.iges "
                     "--out-dir runs/repro_m17306 --allow-open",
        "repro_result": {n: checks["repro"].get(n) for n in ("series60.iges", "kvlcc2_class.iges",
                                                             "KP505_modelScale_250mm.igs")},
        "external_anchor": None,
        "external_anchor_note": "no literature/experimental reference: the two green cases are the "
                                "repository's own write_iges -> read path (self round trip); their only "
                                "measurable external check is the analytic hull volume of "
                                "tensorcad.geomquality.baselines, and the mesh is low by "
                                f"{fid['series60.iges']['rel_deviation_pct']:.2f}% (series60) / "
                                f"{fid['kvlcc2_class.iges']['rel_deviation_pct']:.2f}% (KVLCC2-class), "
                                "a figure the report never states and no gate checks",
        "claim_vs_measured": [
            "doc:71 '4 faces tessellated, 4 refused' vs report KP505 faces_tessellated=0, faces_refused=4",
            "doc:57 table row '4 read, 4 refused' vs report faces_in_file=4, faces_tessellated=0",
            "doc:7 and doc:51 cite runs/gap_t1_iges/evidence.json as the number source: file does not exist",
            "doc:82-84 gold file numbers (164,926 triangles / 30,368 m^3 / 19,568.7 m^3) have no surviving evidence.json",
        ],
        "trust_boundary": "all deliverables are UNTRACKED files in /root/tp-agent/work/m17250; git branch "
                          "agent/m17250 == main 404be46, i.e. no commit in any ref contains them "
                          "(git log --all -- tensorcad/geomquality/iges_trim.py -> empty)",
        "level": 1,
        "level_reason": "anchor function = trimmed IGES/STEP -> watertight mesh; the only third-party "
                        "trimmed file refuses 4/4 faces (0 triangles) and no test drives entity 144/143 "
                        "through import_iges at all, so no executed case of input->mesh exists for the "
                        "trimmed branch; the green G1-G7 belong to the untrimmed patch-set route that "
                        "geomquality.io already provided (io.py:625 read_iges, io.py:652 states type 144 "
                        "out of scope), re-wrapped with the new audit/repair/gate code. Component math is "
                        "verified with closed-form answers (plate/revolution inversion to 1e-12, exact "
                        "areas, box volume 8.0 m^3 to 1e-12) -> L1, not L2.",
        "level_if_scoped_to_untrimmed_route": 2,
    }


def publish() -> None:
    import psycopg

    checks = json.loads((OUT / "checks.json").read_text())
    ev = verdict_evidence()
    sql = ("INSERT INTO tensor_maturity_verdicts (repo, scope, mission_id, level, evidence, verdict_by) "
           "VALUES ('TensorCAD', 'm17250 GAP-CAD-T1-IGES-PIPE', 17306, 1, %s::jsonb, 'board-audit')")
    with psycopg.connect(DSN) as conn, conn.cursor() as cur:
        cur.execute(sql, (json.dumps(ev, ensure_ascii=False),))
        cur.execute("SELECT repo, scope, mission_id, level, verdict_by, created_at "
                    "FROM tensor_maturity_verdicts WHERE mission_id=17306")
        rows = cur.fetchall()
        conn.commit()
    print("inserted:", rows)
    (DOCS / "verdicts.sql").write_text(
        "-- m17306 cross-check verdict for m17250 (anchor GAP-CAD-T1-IGES-PIPE)\n"
        "-- executed against tensorplatform; evidence JSON is the one stored in the table\n"
        "INSERT INTO tensor_maturity_verdicts (repo, scope, mission_id, level, evidence, verdict_by)\n"
        f"VALUES ('TensorCAD', 'm17250 GAP-CAD-T1-IGES-PIPE', 17306, 1, '{json.dumps(ev, ensure_ascii=False)}'::jsonb, 'board-audit');\n")

    markdown = (DOCS / "audit_report.md").read_text()
    png = (OUT / "figure.png").read_bytes()
    payload = {
        "title": "交差核验 m17250 GAP-CAD-T1-IGES-PIPE：trimmed IGES 进口管道实为 0/4 面 → 判级 L1",
        "summary": ("独立复跑 m17250 的 IGES 进口管道：(1) 全部 6 个交付件真实存在但无任何 commit；"
                    "(2) 唯一第三方 trimmed 文件 KP505 4/4 面被拒、0 三角面，报告 doc:71 自述"
                    "「4 faces tessellated」与工件 0 不符；(3) 表 tensor_benchmark_cases 无该锚点基准行（0/46），"
                    "L4 不可能；(4) 两条绿道是既有 untrimmed 路线自往返，体积比解析船体低 1.52%/1.64% 且无门；"
                    "(5) repro 可复现（体积逐位一致）。判级 L1。"),
        "scenario_key": "audit",
        "mission_id": 17306,
        "tags": ["audit", "TensorCAD", "GAP-CAD-T1-IGES-PIPE", "m17250", "maturity-verdict", "L1"],
        "markdown": markdown,
        "figures": [{"name": "figure.png",
                     "caption": "左：IGES 往返后网格体积相对解析船体体积的偏差（红=实测，虚线=±3%）；右：由水密网格反算的 Cb 与基线 Cb。",
                     "png_base64": base64.b64encode(png).decode()}],
    }
    (OUT / "report_payload.json").write_text(json.dumps({k: v for k, v in payload.items() if k != "figures"}, ensure_ascii=False, indent=1))

    def call(path: str, body: dict | None = None, token: str | None = None) -> dict:
        req = urllib.request.Request(API + path, method="POST" if body is not None else "GET")
        req.add_header("Content-Type", "application/json")
        if token:
            req.add_header("Authorization", "Bearer " + token)
        data = json.dumps(body).encode() if body is not None else None
        with urllib.request.urlopen(req, data=data, timeout=180) as r:
            return json.loads(r.read() or b"{}")

    token = call("/auth/login", {"username": ENV["TENSOR_AGENT_USER"],
                                 "password": ENV["TENSOR_AGENT_PASS"]}).get("access_token")
    res = call("/tensor/reports", payload, token)
    print("report:", json.dumps(res)[:400])


if __name__ == "__main__":
    {"measure": measure, "publish": publish}[sys.argv[1]]()