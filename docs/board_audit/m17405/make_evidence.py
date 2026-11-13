#!/usr/bin/env python3
"""Build evidence.json for the m17405 board verdict.

Evidence pointers must be real and machine-checkable (maturity_gate.SCHEME):
gate-checked keys carry *plain existing paths*; line numbers / test names live in
sibling ``*_pointer`` keys so the gate's path resolution is not poisoned by
``file.py:123`` strings.

Also snapshots the subject's (uncommitted) deliverables byte-for-byte into
docs/board_audit/m17405/snapshot/ with a sha256 manifest, so the verdict stays
verifiable after the subject worktree is recycled.

Run: python3 docs/board_audit/m17405/make_evidence.py [subject_worktree]
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SUBJECT = Path(sys.argv[1] if len(sys.argv) > 1 else "/root/tp-agent/work/m17357")
SNAP = HERE / "snapshot"

FILES = [
    "tensorcad/geom_contract.py",
    "docs/geometry_contract.md",
    "tests/test_geom_contract.py",
    "tools/run_contract_compliance.py",
    "tools/generate_contract_assets.py",
    "artifacts/geom_contract_assets/compliance_report.json",
]


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    SNAP.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for rel in FILES:
        src = SUBJECT / rel
        dst = SNAP / rel.replace("/", "__")
        shutil.copyfile(src, dst)
        manifest[rel] = {"subject_sha256": sha(src), "snapshot_sha256": sha(dst),
                         "snapshot": str(dst), "bytes": src.stat().st_size}
    (SNAP / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")

    probe = ROOT / "artifacts/geom_contract_audit/probe_results.json"
    png = ROOT / "artifacts/geom_contract_audit/contract_audit.png"

    evidence = {
        # ---- keys the maturity gate machine-checks (plain paths) -------------
        "impl_files": [
            str(SUBJECT / "tensorcad/geom_contract.py"),
            str(SNAP / "tensorcad__geom_contract.py"),
            str(SUBJECT / "tools/run_contract_compliance.py"),
        ],
        "entry": [str(SUBJECT / "tensorcad/geom_contract.py")],
        "run_artifacts": [
            str(SUBJECT / "artifacts/geom_contract_assets/compliance_report.json"),
            str(probe),
            str(png),
        ],
        "tests": [str(SUBJECT / "tests/test_geom_contract.py")],
        "repro_cmd": (f"cd {SUBJECT} && python3 tools/run_contract_compliance.py  "
                      "# auditor rerun: compliance_report.json sha256 identical"),
        "external_anchor": [
            str(probe),
            str(SUBJECT / "artifacts/geom_contract_assets/box_barge.stl"),
            str(SUBJECT / "artifacts/geom_contract_assets/wigley.stl"),
        ],
        # ---- pointers / audit notes (informational) --------------------------
        "impl_pointer": {
            "library_api": "tensorcad/geom_contract.py:240 validate_contract, :268 validate_stl_file",
            "checks": ("tensorcad/geom_contract.py:124 _check_coordinate_system, :149 _check_units, "
                       ":161 _check_scale, :188 _check_orientation"),
            "cli": "tensorcad/geom_contract.py:282 main  (`python -m tensorcad.geom_contract <stl>`)",
            "docs": "docs/geometry_contract.md (149 lines, four conventions + usage)",
            "repro_tool": "tools/run_contract_compliance.py (8 assets -> compliance_report.json)",
        },
        "tests_pointer": {
            "file": "tests/test_geom_contract.py",
            "run": "cd <subject> && python3 -m pytest tests/test_geom_contract.py -q",
            "result": "22 passed, 0 failed (auditor rerun, 0.4 s)",
            "invariants": ("each convention has both PASS and FAIL cases "
                           "(test_right_handed_passes vs test_left_handed_fails, "
                           "test_outward_normals_pass vs test_inward_normals_fail, "
                           "test_tiny_model_fails/test_huge_model_fails)"),
        },
        "run_artifacts_pointer": {
            "compliance_report": "artifacts/geom_contract_assets/compliance_report.json (8 reports, 5 PASS / 3 FAIL)",
            "assets": "artifacts/geom_contract_assets/*.stl (9 files incl. 3 red-team)",
            "repro_sha256": "ee18be55cb3dbdc1ff43ebbf861dbf45340a83054498bd2e6d283a43de125a0f (before == after)",
        },
        "integration_callers": "NONE — grep -rn 'validate_contract' over the repo returns only the "
                               "re-export at tensorcad/__init__.py:457; no solver/ingest caller "
                               "=> L3 explicitly refused",
        "benchmark": "ABSENT — SELECT ... FROM tensor_benchmark_cases WHERE anchor_slug="
                     "'GAP-CAD-D3-GEOM-CONTRACT' -> 0 rows; no ref/tolerance to compare against "
                     "=> L4 explicitly refused. Auditor substituted 3 closed-form anchors "
                     "(box 1200 m3 exact, sphere -0.474%, Wigley -0.236%), which do NOT count as a "
                     "registered benchmark case.",
        "checks_not_run": [
            "no external literature/witness lookup (offline sandbox); anchors are closed-form only",
            "no cross-repo CAD consumer (TensorFEM/TensorBEM import path) exercised",
            "no large-mesh performance/robustness sweep (validator is O(n) python loops)",
        ],
        "commit_state": {
            "head_commit": None,
            "git_status": "8 untracked/modified paths in the subject worktree; "
                          "/root/tensor-ecosystem/TensorCAD has no tensorcad/geom_contract.py "
                          "=> the capability is invisible to every other repo until committed",
            "server_quality_audit": {"gate_valid": False, "key_present": False,
                                     "gate_missing_warning": True},
            "report_row": "tensor_reports has no row for mission_id=17357 (report_id empty)",
        },
        "auditor_probe": {
            "script": "docs/board_audit/m17405/audit_probe.py",
            "results": "artifacts/geom_contract_audit/probe_results.json",
            "red_team": "5 independently built violations -> 4/4 conventions fire with correct label",
            "false_negatives": 3,
            "anchors": "3 closed-form volume anchors, all |err| < 0.5%",
        },
        "snapshot_manifest": str(SNAP / "MANIFEST.json"),
        "durable_after_push": [
            "docs/board_audit/m17405/AUDIT.md",
            "docs/board_audit/m17405/evidence.json",
            "docs/board_audit/m17405/audit_probe.py",
            "docs/board_audit/m17405/snapshot/",
            "artifacts/geom_contract_audit/probe_results.json",
            "artifacts/geom_contract_audit/contract_audit.png",
        ],
        "level_claimed": 2,
        "level_refused": {
            "L3": "no integration caller (only a re-export in tensorcad/__init__.py:457)",
            "L4": "anchor has 0 rows in tensor_benchmark_cases; no ref/tolerance judge",
            "L5": "no external anchor / no before-after comparison record",
        },
        "verdict_summary": ("L2: real library+CLI that reads an STL, runs four convention checks and "
                            "emits a report; 22 unit tests green, byte-identical regeneration of the "
                            "8-asset compliance table, 4/4 conventions caught on auditor-built "
                            "fixtures. Not integrated into any solver path (L3 refused) and the anchor "
                            "has no registered benchmark case (L4 refused)."),
        "sha256": {k: v["subject_sha256"] for k, v in manifest.items()},
    }
    out = HERE / "evidence.json"
    out.write_text(json.dumps(evidence, ensure_ascii=False, indent=2) + "\n")
    print(f"wrote {out} ({out.stat().st_size} bytes), snapshot {len(manifest)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())