#!/usr/bin/env python3
"""m17402 交差审计落表 + 报告发布：m17350 GAP-CAD-B2-MESH-FORMAT-MATRIX (TensorCAD)。

evidence 全部为本单实测指针（文件:行号 / sha256 / SQL 结果 / 实跑 rc），零推断。
"""
from __future__ import annotations

import base64
import json
import os
import subprocess
import urllib.request
from pathlib import Path

PSQL = ["runuser", "-u", "postgres", "--", "psql", "-d", "tensorplatform", "-tA", "-F", "|", "-c"]
W = Path("/root/tp-agent/work/m17350")
OUT = Path(__file__).resolve().parent


def q(sql: str) -> str:
    return subprocess.run(PSQL + [sql], capture_output=True, text=True).stdout.strip()


def sha(p: Path) -> str:
    import hashlib
    return hashlib.sha256(p.read_bytes()).hexdigest()


bench_rows = q("SELECT count(*) FROM tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-B2-MESH-FORMAT-MATRIX'")
n_all = q("SELECT count(*) FROM tensor_benchmark_cases")
cad_anchor = q("SELECT id||' '||anchor_slug FROM tensor_benchmark_cases WHERE anchor_slug LIKE 'GAP-CAD%'")
head = subprocess.run(["git", "-C", str(W), "rev-parse", "agent/m17350"],
                      capture_output=True, text=True).stdout.strip()
untracked = subprocess.run(["git", "-C", str(W), "status", "--short"],
                           capture_output=True, text=True).stdout.strip()
remote = subprocess.run(["git", "-C", str(W), "ls-remote", "origin"],
                        capture_output=True, text=True).stdout
has_remote = "17350" in remote

evidence = {
    "anchor": "GAP-CAD-B2-MESH-FORMAT-MATRIX",
    "audited_mission": 17350,
    "audited_worktree": f"/root/tp-agent/work/m17350 branch=agent/m17350 head={head[:12]}（=main head，无独立 commit）",
    "impl_files": [
        "tensorcad/meshio_exchange.py:1-182 sha256=f8d8517c92d1efce77df06ca5232e548521baddfc2729adedf05435dee7b72c6："
        "MeshioFormatError(:59)/mesh_to_meshio(:57)/meshio_to_mesh(:80)/write_meshio(:113)/read_meshio(:141)/"
        "round_trip_meshio(:160)；真代码 输入网格→meshio 转换→输出网格，无硬编码结果",
        "tests/test_meshio_format_matrix.py:1-324 sha256=20a0a2d7b1589848c8aaf573b53e23242b02ac2e08d50064170dffb4021a3f4a",
    ],
    "entry": "python -m pytest tests/test_meshio_format_matrix.py -v（repro_cmd.sh:10-11 原文）；"
             "公共 API：from tensorcad.meshio_exchange import round_trip_meshio/write_meshio/read_meshio",
    "run_artifacts": [
        "runs/m17350_meshio_matrix/run_artifacts.json 1043B sha256=29dce299abf2a7b907339718ce22d74a81abd3372a101cfd00aaf5cfda79c088",
        "runs/m17350_meshio_matrix/repro_cmd.sh sha256=b91b6b56c51af732c5ffe0f72c5ced8e7163ea89eb74edce81ca1a0b87c55b7f",
        "artifacts/m17350_meshio_matrix/volume_error_comparison.png 40957B sha256=bdb87c90235e81e769164aa6a7b45c7253c1aced126f52f87b3be49c5144bc61",
        "artifacts/m17350_meshio_matrix/topology_preservation.png 44124B sha256=e149e3eb9e4d649fa627bb1e131310d1f6277838d8f168797c344dd8bd5e4953",
        "⚠️F1：5 个交付文件全部 untracked（git status 全 ??），branch agent/m17350 head=505d10d=main，"
        f"远端分支存在={has_remote}；任务 result 自证 branch=null/head_commit=null、server_quality_audit.gate_valid=false"
        "（TP_AUTO_APPROVE_DEV_ACCEPTANCE 自动放行）→ 幽灵工件，工作树清理即失",
    ],
    "tests": "53/53 PASS（本单实跑 2026-11-13，1.29s，rc=0）：test_round_trip_volume_preserved[12] rel_err<1e-3；"
             "test_round_trip_topology_preserved[12]；test_round_trip_orientation_preserved[12]；"
             "test_exact_topology_for_lossless_formats[9]（atol 1e-12+索引全等）；"
             "test_volume_formats_reject_surface_mesh[6]（pytest.raises MeshioFormatError）；"
             "test_mesh_to_meshio_preserves_data；test_meshio_to_mesh_roundtrip —— 断言均为不变量，非仅执行",
    "integration_callers": "none：grep -rn 'meshio_exchange' --include=*.py（排除模块自身与其测试）0 命中 → L3 不达",
    "benchmark": {
        "registered_case": f"ABSENT —— SELECT count(*) FROM tensor_benchmark_cases WHERE "
                           f"anchor_slug='GAP-CAD-B2-MESH-FORMAT-MATRIX' → {bench_rows} 行（全表 {n_all} 行；"
                           f"唯一 CAD 锚点：{cad_anchor}）；平台裁判无此案号，L4 无法被机器裁定",
        "self_check_recomputed": "本单独立复算（散度定理体积，不复用被测断言）：12/12 支持组合体积相对误差 0.000e+00"
                                "（<1e-15，判据 1e-3）；拓扑 8/8·12/12、513/513·992/992、341/341·678/678 全保；"
                                "OBJ/VTK/Gmsh exact=True；STL 重排但体积/取向保；SU2/CGNS 6/6 正确拒绝 —— 与报告逐位一致",
    },
    "repro_cmd": "cd /root/tp-agent/work/m17350 && python -m pytest tests/test_meshio_format_matrix.py -v —— "
                 "本单实跑 rc=0，53 passed in 1.29s；repro_cmd.sh 的 pip 行（meshio>=5.3）与本机 meshio 5.3.5 相符",
    "external_anchor": "无外部参考值（自指往返保真类任务，无 ref 可对标）；几何锚为解析值 cube V=1.0 / "
                       "sphere 4π/3≈4.18879（实测 −4.121942，法向朝内符号一致）/ wigley 开面半船体（实测 −0.034432）",
    "report": "tensor_reports id=16345 六段（## 1–## 6）正则全命中，figures=2 但 png_base64 为空（F3：图未内嵌，PNG 仅存工作树）",
    "level_rationale": "L2：真代码 输入→计算→输出 + 53 项不变量断言全绿 + 本单独立复算逐位一致。"
                       "不达 L3：无任何集成调用者（grep 0 命中）。不达 L4：(a) tensor_benchmark_cases 无案号（F2）；"
                       "(b) 自指往返无外部参考；(c) 交付物未 commit/未推远端（F1），复现仅在工作树内成立。",
    "verdict_summary": "三件套：run_artifacts ✅存在（但 untracked，F1）；ref vs got ✅自指判据全过且独立复算一致"
                       "（但基准案号缺席，F2）；repro_cmd ✅本单实跑 53/53。报告六段 ✅齐全（图未内嵌 F3）。判级 L2。",
}

sql = (
    "INSERT INTO tensor_maturity_verdicts (repo, scope, mission_id, level, evidence, verdict_by) "
    "VALUES ('TensorCAD', 'm17350 GAP-CAD-B2-MESH-FORMAT-MATRIX', 17402, 2, $json$%s$json$::jsonb, 'board-audit') "
    "ON CONFLICT (repo, scope, verdict_by) DO UPDATE SET level=EXCLUDED.level, evidence=EXCLUDED.evidence, "
    "mission_id=EXCLUDED.mission_id, gate_state='pending', gate_level=NULL, gate_note=NULL "
    "RETURNING id, level, gate_state;" % json.dumps(evidence, ensure_ascii=False)
)
(OUT / "insert_verdict.sql").write_text(sql, encoding="utf-8")
(OUT / "audit_evidence.json").write_text(json.dumps(evidence, ensure_ascii=False, indent=2), encoding="utf-8")
print("INSERT:", subprocess.run(PSQL + [sql], capture_output=True, text=True).stdout.strip())
print("VERIFY:", q("SELECT id,repo,scope,mission_id,level,verdict_by,gate_state FROM tensor_maturity_verdicts "
                   "WHERE repo='TensorCAD' AND scope LIKE 'm17350%' AND verdict_by='board-audit'"))

# ---- 报告发布（六段 + 内嵌被审工件 2 图）----
md = (OUT / "audit_report.md").read_text(encoding="utf-8")
figs = []
for name, cap in [("volume_error_comparison.png", "m17350 工件：12 组合往返体积误差对比（全部 <1e-15，判据线 0.1%）"),
                  ("topology_preservation.png", "m17350 工件：往返前后顶点/面计数对比（全保）")]:
    png = (W / "artifacts/m17350_meshio_matrix" / name).read_bytes()
    figs.append({"name": name, "caption": cap, "png_base64": base64.b64encode(png).decode()})

payload = {
    "title": "m17350 交差审计：TensorCAD meshio 全格式矩阵（GAP-CAD-B2-MESH-FORMAT-MATRIX）—— 判级 L2",
    "summary": "三件套核验：53/53 测试实跑通过、独立复算与报告逐位一致、repro 可执行；但 5 交付文件全部 untracked"
                "（无 commit/无远端分支）且 tensor_benchmark_cases 无本锚点案号 → 判 L2（不达 L3/L4），F1/F2 为返修项。",
    "scenario_key": "GAP-CAD-B2-MESH-FORMAT-MATRIX",
    "mission_id": 17402,
    "tags": ["board-audit", "TensorCAD", "meshio", "maturity-verdict", "L2"],
    "markdown": md,
    "figures": figs,
}
(OUT / "report_payload.json").write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")

tok = json.load(open("/tmp/tok.json"))["token"]
req = urllib.request.Request(os.environ["TENSOR_API_BASE"].rstrip("/") + "/tensor/reports",
                             data=json.dumps(payload).encode(), method="POST")
req.add_header("Content-Type", "application/json")
req.add_header("Authorization", "Bearer " + tok)
with urllib.request.urlopen(req, timeout=120) as r:
    print("REPORT:", r.status, r.read().decode()[:400])
